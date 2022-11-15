from pysam import VariantFile
import argparse


def correct_positions(input_vcf, label_entire_event=False):
    # reads input vcf, adjusts POS, END, TARGET values to account for previous events in the reference
    in_vcf = VariantFile(input_vcf)
    hd = in_vcf.header
    out_vcf = VariantFile(input_vcf[:-4] + '_pos_corrected.vcf', 'w', header=hd)
    # store vcf records in dictionary keyed on start position
    vcf_recs = {}
    for rec in in_vcf.fetch():
        vcf_recs[rec.start] = rec

    # aggregate length of insertions up to a given point in the reference
    # --> to be added to the positions in a given vcf record
    total_ins_len = 0
    for pos in sorted(vcf_recs.keys()):
        # *** shifting the start position of a vcf record automatically shifts the stop position as well
        vcf_recs[pos].start += total_ins_len
        # vcf_recs[pos].stop += total_ins_len
        vcf_recs[pos].info['TARGET'] += total_ins_len
        total_ins_len += vcf_recs[pos].info['SVLEN']

    for rec in vcf_recs.values():
        if not label_entire_event:
            out_vcf.write(rec)
        else:
            # simplified record just reporting start and end positions with respect to entire "A_A*" region
            rec.stop = rec.info['TARGET'] + rec.info['SVLEN']
            rec.info['TARGET'] = -1
            rec.info['SVLEN'] = rec.stop - rec.start
            out_vcf.write(rec)

    out_vcf.close()
    in_vcf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Correct event positions in vcf recording div_dDUPs')
    parser.add_argument('--input_vcf', help='Input vcf from R1 generation step')
    parser.add_argument('--label_entire_event', action='store_true', help='Boolean flag for whether or not'
                                                                          'output vcf should report one record '
                                                                          'per div. repeat with start/end giving '
                                                                          'the overall span of the entire event ('
                                                                          'rather than interval A or a '
                                                                          'specifically)')
    args = parser.parse_args()

    correct_positions(args.input_vcf, args.label_entire_event)
