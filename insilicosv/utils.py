from collections import defaultdict
from contextlib import closing, contextmanager
import dataclasses
from dataclasses import dataclass
from enum import Enum
import logging
import os
import os.path
import random
from typing_extensions import Optional, override
from intervaltree import IntervalTree
import pysam
from sortedcontainers import SortedSet  # type: ignore
from copy import deepcopy

logger = logging.getLogger(__name__)


def if_not_none(a, b):
    return a if a is not None else b


def chk(cond, msg="Error", error_type='runtime'):
    errors = {'runtime': RuntimeError, 'value': ValueError, 'type': TypeError,
              'index': IndexError, 'syntax': SyntaxError, 'file not found': FileNotFoundError}
    if not cond:
        raise errors[error_type](f"insilicoSV error: {msg}")


def generate_seq(length):
    base_map = {1: "A", 2: "T", 3: "G", 4: "C"}
    return ''.join([base_map[random.randint(1, 4)] for _ in range(length)])


def complement(seq):
    output = ""
    base_complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    for base in seq.upper():
        if base in base_complements:
            output += base_complements[base]
        else:
            output += base

    return output


def reverse_complement(seq):
    return complement(seq[::-1])


def remove_file(path):
    if os.path.exists(path):
        os.remove(path)


def divergence(seq, divergence_prob):
    # apply random base flips to input sequence
    p = float(divergence_prob)
    assert 0 <= p <= 1
    if p == 0:
        return seq

    def mutate_base(b):
        if b not in 'TCGA':
            return b
        return random.choice([allele for allele in 'TCGA' if allele != b])

    return ''.join([b if random.random() > p else mutate_base(b) for b in seq.upper()])


def is_readable_file(fname):
    try:
        with open(fname) as test_f:
            pass
        return True
    except:
        return False


def as_list(val):
    if isinstance(val, str):
        val = [val]
    return val


def is_valid_closed_int_range(val):
    return (
            isinstance(val, list) and
            len(val) == 2 and
            isinstance(val[0], int) and
            isinstance(val[1], int) and
            0 <= val[0] <= val[1])


def is_list_of(type_, val):
    return isinstance(val, list) and all(isinstance(v, type_) for v in val)


@dataclass(order=True, frozen=True)
class Locus:
    chrom: str
    pos: int

    def shifted(self, delta: int):
        return Locus(chrom=self.chrom, pos=self.pos + delta)


@dataclass(order=True, frozen=True)
class Region:
    """A specific, contiguous genomic region: (chrom, start, end).

    The region includes `start` but not `end`.  An empty region with start==end represents
    an insertion point before `start`.  Insertions at the same point are ordered by `order_key`.
    The coordinates `start` and `end` are 0-based.
    """

    chrom: str
    start: int
    end: int
    order_key: tuple = ()

    # Optional extra info that can be associated with Regions representing
    # ROIs or pieces of ROIs:

    # region type (e.g. L1, Alu, etc)
    kind: str = ''

    # additional per-region data, such as repeat unit info for tandem repeats
    data: str = 0
    motif: str = ''

    # if this region is derived from an ROI, start/end of the original ROI
    orig_start: int = -1
    orig_end: int = -1

    # keep the operation in the region for the source regions of overlapping SVs
    sv: Optional['SV'] = None

    def __post_init__(self):
        assert self.chrom
        assert self.start <= self.end

    def is_empty(self):
        return self.start == self.end

    def length(self):
        return self.end - self.start

    def replace(self, **kw):
        """Return a copy of self with fields replaced according to `kw`"""
        return dataclasses.replace(self, **kw)

    def shifted(self, delta: int):
        return self.replace(start=self.start + delta, end=self.end + delta)

    def padded(self, padding: int):
        return self.replace(start=self.start - padding, end=self.end + padding)


# end: class Region

@dataclass(frozen=True)
class RegionFilter:
    region_kinds: Optional[tuple[str, ...]] = None
    region_length_range: tuple[Optional[int], Optional[int]] = (None, None)

    def satisfied_for(self, region) -> bool:
        if (self.region_kinds is not None and
                (not region.kind or
                 not any((region_kinds.upper() == 'ALL' and region.kind != '_reference_') or
                         region_kinds in region.kind
                         for region_kinds in self.region_kinds))):
            return False
        return True


@dataclass(frozen=True)
class TandemRepeatRegionFilter(RegionFilter):
    min_num_repeats: int = 0

    @override
    def satisfied_for(self, region: Region) -> bool:
        if not super().satisfied_for(region):
            return False
        try:
            repeat_unit_length: int = int(region.data)
        except ValueError:
            return False
        return repeat_unit_length * self.min_num_repeats <= region.length()


class OverlapMode(Enum):
    """Types of spatial relationship between an anchor (SV sub-region) and an ROI."""

    PARTIAL = "partial"  # exactly one endpoint of anchor in ROI
    CONTAINED = "contained"  # both endpoints of anchor in ROI
    EXACT = "exact"  # anchor exactly matches the ROI
    CONTAINING = "containing"  # both endpoints of ROI strictly inside the anchor
    TERMINAL = "terminal"  # One breakpoint at a chromosome extremity
    CHROM = "whole-chromosome"  # The breakpoints are a chromosome extremities (only for DUP and DEL), for DUP creates chromosome copies


class RegionSet:
    """A collection of genomic regions"""

    chrom2itree: dict[str, IntervalTree]
    num_hap_tree: int

    def __init__(self, regions=None, allow_hap_overlap=False):
        regions = regions or []
        self.num_hap_tree = 3 if allow_hap_overlap else 1
        chrom2regions = defaultdict(list)
        for region in regions:
            chrom2regions[region.chrom].append(region)

        # One tree for each haplotype and one for the combination to check homozygous variants
        self.chrom2itree = defaultdict(lambda: defaultdict(IntervalTree))
        for chrom, chrom2region in chrom2regions.items():
            if len(chrom2region) > 100000:
                logger.debug(f'RegionSet init: {chrom=} {len(chrom2region)=}')
            for hap in range(self.num_hap_tree):
                # Padding to allow insertion target regions
                self.chrom2itree[chrom][hap] = IntervalTree.from_tuples((region.start - 0.1, region.end + 0.1, region)
                                                                        for region in chrom2region)

    def __contains__(self, region):
        overlap = []
        for hap in range(self.num_hap_tree):
            overlap += list(self.chrom2itree[region.chrom].overlap(region.start, region.end))
        return any([(interval.begin == region.start and interval.end == region.end) for interval in overlap])

    def strictly_contains_point(self, point, chrom):
        overlap = list(self.chrom2itree[chrom][0].at(point))
        for interval in overlap:
            if interval.data.start < point < interval.data.end:
                return True
        return False

    @staticmethod
    def from_beds(bed_paths, to_region_set, verbose=False):
        regions = []
        for bed_path in bed_paths:
            logger.info(f'Reading bed file {bed_path}')
            with open(bed_path) as bed:
                for line_num, line in enumerate(bed):
                    if (verbose and (line_num % 500000) == 0):
                        logger.debug(f'line {line_num}')
                    loc = f'file {bed_path}, line {line_num}'
                    if line.startswith('#') or line.isspace():
                        continue
                    fields = line.strip().split()
                    chk(len(fields) >= 3,
                        f'{loc}: too few fields in line in the BED file {bed_path}', error_type='value')
                    chrom, start_str, end_str = fields[:3]
                    kind = fields[3] if len(fields) > 3 else 'NA'
                    chk(all((chrom, start_str, end_str)),
                        f'{loc}: empty value in first three columns in the BED file {bed_path}', error_type='value')
                    try:
                        start, end = int(start_str), int(end_str)
                    except ValueError:
                        chk(False, f'{loc}: invalid start or end of region in the BED file {bed_path}',
                            error_type='value')
                    chk(start < end, f'{loc}: region start must be less than end in the BED file {bed_path}',
                        error_type='value')
                    chk(0 <= start, f'{loc}: region start must be non-negative in the BED file {bed_path}',
                        error_type='value')
                    motif = fields[4] if len(fields) >= 5 else ''
                    data = len(motif)
                    regions.append(Region(chrom=chrom, start=start,
                                          end=end, kind=kind, data=data,
                                          motif=motif, orig_start=start, orig_end=end))

        region_set = regions
        if to_region_set:
            logger.info(f'Constructing Interval Tree from {len(regions)} regions...')
            region_set = RegionSet(regions)
            logger.info(f'Constructed Interval Tree from {len(regions)} regions.')
        return region_set

    @staticmethod
    def from_vcf(vcf_path):
        regions = []
        with closing(pysam.VariantFile(vcf_path)) as vcf_file:
            for vcf_rec in vcf_file.fetch():
                vcf_info = dict(vcf_rec.info)
                kind = 'DEFAULT'
                if 'REGION_TYPE' in vcf_info:
                    kind = vcf_info['REGION_TYPE']
                regions.append(Region(chrom=vcf_rec.chrom, start=vcf_rec.start, end=vcf_rec.stop,
                                      orig_start=vcf_rec.start, orig_end=vcf_rec.stop, kind=kind))
                if 'TARGET' in vcf_info and isinstance(vcf_info['TARGET'], int):
                    target_chrom = vcf_info.get('TARGET_CHROM', vcf_rec.chrom)
                    target_start = vcf_info['TARGET'] - 1
                    regions.append(Region(chrom=target_chrom, start=target_start, end=target_start,
                                          orig_start=target_start, orig_end=target_start, kind=kind))

        return RegionSet(regions)

    @staticmethod
    def from_fasta(fasta_path, filter_small_chr, region_kind, allow_hap_overlap):
        with pysam.FastaFile(fasta_path) as fasta_file:
            regions = []
            for chrom, chrom_length in zip(fasta_file.references, fasta_file.lengths):
                if chrom_length < filter_small_chr: continue
                regions.append(Region(chrom=chrom, start=0, end=chrom_length,
                                      kind=region_kind,
                                      orig_start=0, orig_end=chrom_length))
            return RegionSet(regions, allow_hap_overlap=allow_hap_overlap)

    def get_region_list(self, hap=0):
        return [ival.data for chrom_itree in self.chrom2itree.values() for ival in chrom_itree[hap]]

    def filtered(self, region_filter):
        """Construct a RegionSet containing regions from self that meet given filter"""

        # TODO: use numpy to filter regions faster

        def satisfies_filter(region):
            return region_filter.satisfied_for(region)

        # Blacklist affect the three haplotypes
        return RegionSet(filter(satisfies_filter, self.get_region_list()))

    def add_region_set(self, other_region_set):
        for chrom, other_chrom_itree_list in other_region_set.chrom2itree.items():
            for hap, other_chrom_itree in other_chrom_itree_list.items():
                self.chrom2itree[chrom][hap].update(other_chrom_itree)

    def add_region(self, region, sv=None, allow_hap_overlap=None):
        aux_region = deepcopy(region)
        if sv:
            aux_region = aux_region.replace(sv=sv)
        self.add_region_set(RegionSet([aux_region], allow_hap_overlap=allow_hap_overlap))

    def chop(self, sv_region, genotype):
        # Remove the parts of intervals overlapping sv_region.
        for hap in range(self.num_hap_tree):
            # Only chop the intervals of the tree on the same haplotype
            if hap != self.num_hap_tree - 1 and not genotype[hap]: continue

            chrom_itree = self.chrom2itree[sv_region.chrom][hap]

            def adjust_region(ival, is_begin):
                if is_begin:
                    new_region = ival.data.replace(end=round(sv_region.start))
                else:
                    new_region = ival.data.replace(start=round(sv_region.end))
                if new_region.start == new_region.end:
                    return None
                return new_region

            # Pad to fully remove insertion target
            chrom_itree.chop(sv_region.start - 0.1, sv_region.end + 0.1, datafunc=adjust_region)


# end class RegionSet

def percent_N(seq):
    return 0 if len(seq) == 0 else (seq.count('N') + seq.count('n')) / len(seq)


def pairwise(iterable):
    # pairwise('ABCD') → [('A', 'B'), ('B', 'C'), ('C', 'D')]
    iterator = iter(iterable)
    a = next(iterator, None)
    for b in iterator:
        yield a, b
        a = b


def has_duplicates(it):
    items = set()
    for item in it:
        if item in items:
            return True
        items.add(item)
    return False


@contextmanager
def error_context(*args):
    context_messages = tuple(map(str, args))
    try:
        yield
    except RuntimeError as e:
        e.args += context_messages
        raise