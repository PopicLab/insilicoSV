from collections import defaultdict, Counter, namedtuple
from contextlib import closing
import logging
import os
import os.path
import re
from pysam import VariantFile

from insilicosv import utils
from insilicosv.utils import Region, Locus
from insilicosv.sv_defs import Operation, Transform, TransformType, BreakendRegion
from insilicosv.variant_set import get_vcf_header_infos

logger = logging.getLogger(__name__)

class StatsCollector:

    def __init__(self, chroms, chrom_lengths):
        assert len(chroms) == len(chrom_lengths)
        self.num_heterozygous = 0
        self.num_homozygous = 0
        self.placed_svs = 0
        self.total_svs = 0
        self.chroms = chroms
        self.chrom_lengths = chrom_lengths
        self.len_frags_chr = defaultdict(int)  # Lengths of altered fragments within chromosome
        self.sv_types = Counter()
        self.region_types = Counter()

    def get_info(self, svs):
        """
        collects all information for stats file after all edits are completed
        """
        self.total_svs = len(svs)
        
        self.min_region_len = None
        self.max_region_len = None

        tot_frags_len: int = 0
        tot_frags_count = 0

        for sv in svs:
            if not sv.is_placed():
                continue
            self.placed_svs += 1
            sv_type = sv.info.get('SVTYPE', 'UNKNOWN')
            self.sv_types[sv_type] += 1
            if (sv.roi is not None) and (sv.roi.kind != '_reference_'):
                self.region_types[sv.roi.kind] += 1
            assert sv.genotype is not None
            if sv.genotype[0] and sv.genotype[1]:
                self.num_homozygous += 1
            else:
                self.num_heterozygous += 1

            for src_reg in set([operation.source_region for operation in sv.operations
                                if operation.source_region is not None]):
                self.len_frags_chr[src_reg.chrom] += src_reg.length()

                tot_frags_len += src_reg.length()
                tot_frags_count += 1

                if self.min_region_len is None or src_reg.length() < self.min_region_len:
                    self.min_region_len = src_reg.length()
                if self.max_region_len is None or src_reg.length() > self.max_region_len:
                    self.max_region_len = src_reg.length()
        self.avg_len = tot_frags_len // tot_frags_count if tot_frags_count != 0 else 0

    def export_data(self, fileout):
        """
        Exports all collected data to entered file
        fileout: Location to export stats file
        """
        def write_item(fout, name, item, prefix=""):
            fout.write("{}{}: {}\n".format(prefix, str(name), str(item)))

        os.makedirs(os.path.dirname(fileout), exist_ok=True)
        with open(fileout, "w") as fout:
            fout.write("===== Overview =====\n")
            write_item(fout, "SVs successfully simulated", str(self.placed_svs) + "/" + str(self.total_svs))
            for sv_type in self.sv_types:
                write_item(fout, sv_type, self.sv_types[sv_type], prefix="\t- ")
            write_item(fout, "Homozygous SVs", self.num_homozygous)
            write_item(fout, "Heterozygous SVs", self.num_heterozygous)
            write_item(fout, "Average length of impacted reference regions", self.avg_len)
            write_item(fout, "Min length of impacted reference region", self.min_region_len)
            write_item(fout, "Max length of impacted reference region", self.max_region_len)
            if self.region_types:
                write_item(fout, "NUmber of SVs per ROI types", self.region_types)
            for chrom, chrom_length in zip(self.chroms, self.chrom_lengths):
                fout.write(f"\n===== {chrom} =====\n")
                write_item(fout, "Length of sequence", chrom_length)
                write_item(fout, "Total impacted length of reference chromosome",
                           self.len_frags_chr[chrom])
# end: class StatsCollector

PafRecord = namedtuple('PafRecord', [
    "query_name",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "target_name",
    "target_length",
    "target_start",
    "target_end",
    "residue_matches",
    "alignment_block_length",
    "mapping_quality",
    "tags",
])

class OutputWriter:

    def __init__(self, svs, reference, chrom_lengths, output_path, config):
        self.svs = svs
        self.reference = reference
        self.chrom_lengths = chrom_lengths
        self.output_path = output_path
        self.config = config

    def output_haps(self):
        if self.config.get('output_no_haps', False):
            logger.warning('Skipping haps output')
            return
        for hap_index, hap_fa in enumerate(['sim.hapA.fa', 'sim.hapB.fa']):
            self.output_hap(os.path.join(self.output_path, hap_fa), hap_index)

    def output_hap(self, hap_fa, hap_index):
        paf_records = []
        hap_chrom_lengths = {}
        with open(hap_fa, 'w') as sim_fa:
            chrom2operations: dict[str, list[Operation]] = self.get_chrom2operations(hap_index)
            for chrom, chrom_length in zip(self.reference.references, self.reference.lengths):
                hap_str = ['hapA', 'hapB'][hap_index]
                hap_chrom = f'{chrom}_{hap_str}'
                sim_fa.write(f'>{hap_chrom}\n')

                hap_chrom_pos = 0
                for operation in chrom2operations[chrom]:
                    # paf format: https://cran.r-project.org/web/packages/pafr/vignettes/Introduction_to_pafr.html
                    paf_mapq = 60
                    for _ in range(operation.transform.n_copies):
                        if operation.transform_type == TransformType.DEL:
                            assert operation.source_region is not None
                            paf_rec = PafRecord(
                                query_name=hap_chrom, query_length=None,
                                query_start=hap_chrom_pos, query_end=hap_chrom_pos,
                                strand='+',
                                target_name=operation.source_region.chrom,
                                target_length=self.chrom_lengths[operation.source_region.chrom],
                                target_start=operation.source_region.start,
                                target_end=operation.source_region.end,
                                residue_matches=0, alignment_block_length=operation.source_region.length(),
                                mapping_quality=paf_mapq,
                                tags=f'cg:Z:{operation.source_region.length()}D')
                        else:
                            if operation.source_region is not None:
                                seq = self.reference.fetch(
                                    reference=operation.source_region.chrom,
                                    start=operation.source_region.start,
                                    end=operation.source_region.end)
                                paf_rec = PafRecord(
                                    query_name=hap_chrom, query_length=None,
                                    query_start=hap_chrom_pos, query_end=hap_chrom_pos + len(seq),
                                    strand='+',
                                    target_name=operation.source_region.chrom,
                                    target_length=self.chrom_lengths[operation.source_region.chrom],
                                    target_start=operation.source_region.start,
                                    target_end=operation.source_region.end,
                                    residue_matches=len(seq), alignment_block_length=len(seq),
                                    mapping_quality=paf_mapq,
                                    tags=f'cg:Z:{len(seq)}M')
                            else:
                                assert operation.novel_insertion_seq is not None
                                seq = operation.novel_insertion_seq
                                paf_rec = PafRecord(
                                    query_name=hap_chrom, query_length=None,
                                    query_start=hap_chrom_pos, query_end=hap_chrom_pos + len(seq),
                                    strand='+',
                                    target_name=operation.target_region.chrom,
                                    target_length=self.chrom_lengths[operation.target_region.chrom],
                                    target_start=operation.target_region.start,
                                    target_end=operation.target_region.end,
                                    residue_matches=0, alignment_block_length=len(seq),
                                    mapping_quality=paf_mapq,
                                    tags=f'cg:Z:{len(seq)}I')

                            if operation.transform_type == TransformType.INV:
                                seq = utils.reverse_complement(seq)
                                paf_rec = paf_rec._replace(strand='-')
                            if (operation.transform.divergence_prob > 0 or
                                operation.transform.replacement_seq is not None):
                                seq_orig = seq
                                seq = utils.if_not_none(
                                    operation.transform.replacement_seq,
                                    utils.divergence(seq, operation.transform.divergence_prob))
                                assert len(seq) == len(seq_orig)
                            sim_fa.write(seq)
                            hap_chrom_pos += len(seq)

                        # end: if operation.transform_type == TransformType.DEL:

                        assert paf_rec is not None

                        if operation.op_info and operation.op_info.get('sv_id'):
                            sv_id = operation.op_info['sv_id']
                            paf_rec = paf_rec._replace(tags=paf_rec.tags + f'\tsv:Z:{sv_id}')
                        paf_records.append(paf_rec)
                    # end: for _ in range(operation.transform.n_copies):
                # end: for operation in chrom2operations[chrom]

                hap_chrom_lengths[hap_chrom] = hap_chrom_pos
                sim_fa.write('\n')
            # end: for chrom, chrom_length in zip(...)
        # end: with open(hap_fa, 'w') as sim_fa
        if self.config.get('output_paf', False):
            hap_paf = hap_fa.replace('.fa', '.paf')

            new_paf_records = []
            for paf_rec_num, paf_rec in enumerate(paf_records):
                if (paf_rec.target_start == paf_rec.target_end and
                    paf_rec.query_start == paf_rec.query_end):
                    continue
                new_paf_records.append(paf_rec)

            with open(hap_paf, 'w') as sim_paf:
                for paf_rec in new_paf_records:
                    paf_rec = paf_rec._replace(query_length=hap_chrom_lengths[paf_rec.query_name])
                    sim_paf.write('\t'.join(map(str, paf_rec)) + '\n')

    # end: def output_hap(self, hap_fa, hap_index)

    def get_chrom2operations(self, hap_index):
        """For each chromosome, make a list of operations targeting that chromosome,
        sorted along the chromosome.  Add identity operations for regions not touched
        by SVs."""
        chrom2operations = defaultdict(list)
        target_region2transform = {}
        for sv in self.svs:
            assert sv.genotype is not None
            if sv.genotype[hap_index]:
                for operation in sv.operations:
                    assert operation.target_region is not None
                    if operation.is_in_place:
                        if operation.target_region in target_region2transform:
                            assert (operation.transform ==
                                    target_region2transform[operation.target_region])
                            continue
                        target_region2transform[operation.target_region] = operation.transform

                    if operation.op_info is None:
                        operation.op_info = dict()
                    operation.op_info['sv_id'] = sv.sv_id
                    chrom2operations[operation.target_region.chrom].append(operation)

        for chrom, chrom_length in self.chrom_lengths.items():
            target_regions = sorted([Region(chrom=chrom, start=0, end=0)] + 
                                    [operation.target_region for operation in chrom2operations[chrom]] +
                                    [Region(chrom=chrom, start=chrom_length, end=chrom_length)])

            intersv_ops = [Operation(
                transform=Transform(transform_type=TransformType.IDENTITY,
                                    is_in_place=True),
                source_breakend_region=BreakendRegion(0,1),
                placement=[Locus(chrom, target_region1.end),
                           Locus(chrom, target_region2.start)])
                for target_region1, target_region2 in utils.pairwise(target_regions)
                if target_region2.start > target_region1.end]
            
            chrom2operations[chrom].extend(intersv_ops)
            chrom2operations[chrom].sort(key=lambda operation: operation.target_region)
        return chrom2operations

    def output_novel_insertions(self):
        with open(os.path.join(self.output_path, 'sim.novel_insertions.fa'), 'w') as (
                sim_novel_insertions_fa):
            for sv in self.svs:
                for op_num, operation in enumerate(sv.operations):
                    if operation.novel_insertion_seq:
                        seq_id = (
                            f'{sv.sv_id}_{op_num}_{operation.target_region.chrom}_'
                            f'{operation.target_region.start}')
                        sim_novel_insertions_fa.write(f'>{seq_id}\n')
                        sim_novel_insertions_fa.write(f'{operation.novel_insertion_seq}\n')

    def output_vcf(self):
        vcf_path = os.path.join(self.output_path, 'sim.vcf')
        with open(vcf_path, "w") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            for chrm, chrm_len in zip(self.reference.references,
                                      self.reference.lengths):
                vcf.write("##contig=<ID=%s,length=%d>\n" % (chrm, chrm_len))
            vcf.write("#%s\n" % "\t".join(["CHROM", "POS", "ID",
                                           "REF", "ALT", "QUAL", "FILTER", "INFO",
                                           "FORMAT", "SAMPLE"]))

        with closing(VariantFile(vcf_path)) as vcf_file:

            for vcf_header_info in get_vcf_header_infos():
                vcf_file.header.info.add(**vcf_header_info)

            vcf_file.header.formats.add('GT', number=1, type='String', description="Genotype")

            with closing(VariantFile(vcf_path, 'w', header=vcf_file.header)) as vcf_out_file:

                vcf_records: list[dict] = []

                for sv in self.svs:
                    vcf_records.extend(sv.to_vcf_records(self.config))
                assert not utils.has_duplicates(vcf_rec['id'] for vcf_rec in vcf_records)
                for vcf_rec in sorted(vcf_records, key=lambda rec: (rec['contig'], rec['start'])):
                    self.check_vcf_record(vcf_rec)
                    vcf_record = vcf_out_file.header.new_record(**vcf_rec)  # type: ignore
                    assert len(vcf_record.samples) == 1
                    vcf_record.samples[0].phased = True
                    vcf_out_file.write(vcf_record)
            # end: with closing(VariantFile(vcf_path, 'w', header=vcf_file.header)) as vcf_out_file
        # end: with closing(VariantFile(vcf_path)) as vcf_file
    # end: def output_vcf(self)

    def check_vcf_record(self, vcf_rec):
        assert set(vcf_rec.keys()) >= {'id', 'contig', 'start', 'stop', 'qual',
                                       'filter', 'alleles', 'samples'}
        assert not any(c.isspace() or c == ';' for c in vcf_rec['id'])
        info = vcf_rec.get('info')
        assert isinstance(info, (type(None), dict))
        for info_key, info_val in vcf_rec.get('info', {}).items():
            if isinstance(info_val, (str, list)):
                for val in utils.as_list(info_val):
                    assert re.search(r'[\s;,=]', val) is None

    def output_svops_bed(self):
        if not self.config.get('output_svops_bed', False):
            return

        with open(os.path.join(self.output_path, 'sim.svops.bed'), 'w') as out:
            for sv in self.svs:
                for operation in sv.operations:
                    op_str = '_'.join([sv.sv_id, sv.info['SVTYPE'], operation.transform_type.value])
                    if operation.op_info and 'symbol' in operation.op_info:
                        op_str += ('_' + operation.op_info['symbol'].name)

                    out.write('\t'.join(map(str,
                                            [operation.target_region.chrom,
                                             operation.target_region.start,
                                             max(operation.target_region.end,
                                                 operation.target_region.start+1),
                                             op_str]))
                              + '\n')

                    if (operation.source_region is not None and
                        operation.source_region != operation.target_region):
                        out.write('\t'.join(map(str,
                                                [operation.source_region.chrom,
                                                 operation.source_region.start,
                                                 max(operation.source_region.end,
                                                     operation.source_region.start),
                                                 op_str + '_src']))
                                  + '\n')

    def output_stats(self):
        stats_collector = StatsCollector(chroms=list(self.reference.references),
                                         chrom_lengths=list(self.reference.lengths))
        stats_collector.get_info(self.svs)
        stats_collector.export_data(os.path.join(self.output_path, 'sim.stats.txt'))

# end: class OutputWriter
