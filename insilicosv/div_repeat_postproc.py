import argparse
from pysam import VariantFile
from collections import defaultdict


def correct_positions_div(input_vcf):
    # reads input vcf, adjusts POS, END, TARGET values to account for previous events in the reference
    in_vcf = VariantFile(input_vcf)
    hd = in_vcf.header
    file_suffix = '_decoy_divrepeat_intervals.vcf'
    out_vcf = VariantFile(input_vcf[:-4] + file_suffix, 'w', header=hd)
    # store vcf records in dictionary keyed on start position
    vcf_recs = defaultdict(dict)
    for rec in in_vcf.fetch():
        vcf_recs[rec.chrom][rec.start] = rec

    # aggregate length of insertions up to a given point in the reference (to be added to the positions in a given vcf record)
    for chrom in vcf_recs.keys():
        total_ins_len = 0
        for pos in sorted(vcf_recs[chrom].keys()):
            evt = vcf_recs[chrom][pos]
            if evt.id == 'div_dDUP':
                evt.start += total_ins_len
                evt.info['TARGET'] += total_ins_len
                total_ins_len += evt.info['SVLEN']
        for rec in vcf_recs[chrom].values():
            # output separate records for each A and A* ("B") to serve as avoid intervals for further simulation
            # ** assuming vcf will only contain SVs of type div_dDUP
            rec_A = rec.copy()
            rec_B = rec.copy()
            rec_A.id = 'div_dDUP_A'
            rec_A.alts = ('div_dDUP_A',)
            rec_A.info['SVTYPE'] = 'div_dDUP_A'
            rec_A.info['TARGET'] = -1
            # setting the above unsets rec.stop, so need to reset it to the original
            rec_A.stop = rec.stop
            # if the event is flipped, then need to shift the _A interval by one interval length
            # to compensate for the insertion at the target locus, that will in this case happen before the source
            if rec.info['TARGET'] < rec.start:
                rec_A.start += rec.info['SVLEN']
            rec_B.id = 'div_dDUP_B'
            rec_B.alts = ('div_dDUP_B',)
            rec_B.start = rec_B.info['TARGET']
            rec_B.stop = rec_B.info['TARGET'] + rec.info['SVLEN']
            rec_B.info['SVTYPE'] = 'div_dDUP_B'
            rec_B.info['TARGET'] = -1
            out_vcf.write(rec_A)
            out_vcf.write(rec_B)

    out_vcf.close()
    in_vcf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Utility script for multi-staged simulation procedure of divergent repeats')
    parser.add_argument('--div_dDUP_vcf', help='Input vcf from R1 generation step')
    args = parser.parse_args()
    correct_positions_div(args.div_dDUP_vcf)
