#!/usr/bin/env python3

"""Make plots summarizing the variants from an insilicoSV simulation.
"""

import argparse
from collections import defaultdict
from contextlib import closing
from dataclasses import dataclass
import os
import os.path
import logging
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import pysam

from plotnine import (
    ggplot, geom_point, aes, facet_wrap, facet_grid,
    ggtitle, theme, scale_y_log10, save_as_pdf_pages)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


@dataclass(order=True, frozen=True)
class Region:
    chrom: str
    start: int
    end: int

    def __post_init__(self):
        assert 0 <= self.start <= self.end

    def length(self) -> int:
        return self.end - self.start

def pairwise(iterable):
    # pairwise('ABCD') → [('A', 'B'), ('B', 'C'), ('C', 'D')]
    iterator = iter(iterable)
    a = next(iterator, None)
    for b in iterator:
        yield a, b
        a = b

def run_plot_sim_stats(args):

    svid2source_regions: dict[str, list[Region]] = defaultdict(list)

    for hap in ('hapA', 'hapB'):
        with open(f'{args.insilicosv_out_dir}/sim.{hap}.paf') as sim_paf:
            for line in sim_paf:
                fields = line.strip().split('\t')
                for tag in fields[12:]:
                    tag_name, tag_type, tag_val = tag.split(':')
                    if tag_name == 'sv':
                        svid2source_regions[tag_val].append(Region(chrom=fields[5],
                                                                   start=int(fields[7]),
                                                                   end=int(fields[8])))

    plots = []
    
    part_records = []

    vcf_fname = f'{args.insilicosv_out_dir}/sim.vcf'

    with closing(pysam.VariantFile(vcf_fname)) as vcf_file:
        for rec in vcf_file.fetch():
            svid = rec.id
            svtype = rec.info['SVTYPE']
            if 'PARENT_SVID' in rec.info:
                svid = rec.info['PARENT_SVID']
                svtype = rec.info['PARENT_SVTYPE']

            if 'TARGET' in rec.info:
                svid2source_regions[svid].append(Region(chrom=rec.info['TARGET_CHROM'],
                                                        start=rec.info['TARGET']-1,
                                                        end=rec.info['TARGET']-1))
            
    with closing(pysam.VariantFile(vcf_fname)) as vcf_file:
        svids_seen = set()
        for rec in vcf_file.fetch():
            svid = rec.id
            svtype = rec.info['SVTYPE']
            if 'PARENT_SVID' in rec.info:
                svid = rec.info['PARENT_SVID']
                svtype = rec.info['PARENT_SVTYPE']

            if svid in svids_seen:
                continue
            svids_seen.add(svid)
            if not args.svtypes or svtype in args.svtypes:

                source_regions = sorted(set(svid2source_regions[svid]))
                print(f'{source_regions=}')
                for source_region in source_regions:
                    if source_region.length() > 0:
                        part_records.append(dict(svtype=svtype,
                                                 part_type='component',
                                                 part_len=source_region.length(),
                                                 chrom=source_region.chrom,
                                                 vset=rec.info['VSET']))

                for region1, region2 in pairwise(source_regions):
                    if region1.chrom == region2.chrom:
                        disp_region = Region(chrom=region1.chrom, start=region1.end, end=region2.start)
                        if disp_region.length() > 0 and disp_region not in source_regions:
                            part_records.append(dict(svtype=svtype,
                                                     part_type='dispersion',
                                                     part_len=disp_region.length(),
                                                     chrom=disp_region.chrom,
                                                     vset=rec.info['VSET']))

    if not part_records:
        raise RuntimeError('Nothing to plot!')

    parts_df = pd.DataFrame(part_records)

    if not args.parts_tsv:
        args.parts_tsv = vcf_fname.replace('.vcf', '.parts.tsv')

    if args.parts_tsv:
        parts_df.to_csv(args.parts_tsv, sep='\t', index=False)

    for vset, vset_df in parts_df.groupby('vset'):
        for part_type, part_df in vset_df.groupby('part_type'):
            plots.append(
                ggplot(data=part_df.rename(columns={'part_len': f'{part_type}_len'}))
                + ggtitle(f'vset{vset} {part_type} lengths')
                + theme(figure_size=(args.figsize[0], args.figsize[1]))
                + geom_point(mapping=aes(x="svtype", y=f"{part_type}_len"), position="jitter", alpha=0.2)
                + scale_y_log10()
            )

    if not args.report_pdf:
        args.report_pdf = vcf_fname.replace('.vcf', '.report.pdf')
    save_as_pdf_pages(plots=plots,
                      filename=args.report_pdf)
    logger.info(f'Saved report to {args.report_pdf}')

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--insilicosv_out_dir', required=True, help='insilicoSV output directory')
    parser.add_argument('--svtypes', nargs='+', default=[], help='SV types to include in the plot')

    parser.add_argument('--figsize', type=int, nargs=2, default=[8, 6], help='figure size for the report')

    parser.add_argument('--parts_tsv', help='output file for variant parts')
    parser.add_argument('--report_pdf', help='output file for the report (pdf)')

    return parser.parse_args()
    
if __name__ == '__main__':
    run_plot_sim_stats(parse_args())
