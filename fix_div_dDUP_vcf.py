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
        evt = vcf_recs[pos]
        # *** shifting the start position of a vcf record automatically shifts the stop position as well
        evt.start += total_ins_len
        # in case of mixed vcf (div_dDUPs and other types) only want to shift events by cumulative div_dDUP length
        if evt.id == 'div_dDUP':
            evt.info['TARGET'] += total_ins_len
            total_ins_len += evt.info['SVLEN']

    for rec in vcf_recs.values():
        if not label_entire_event:
            out_vcf.write(rec)
        elif rec.id == 'div_dDUP':
            # simplified record just reporting start and end positions with respect to entire "A_A*" region
            rec.stop = rec.info['TARGET'] + rec.info['SVLEN']
            rec.info['TARGET'] = -1
            rec.info['SVLEN'] = rec.stop - rec.start
            out_vcf.write(rec)

    out_vcf.close()
    in_vcf.close()


def incorporate_events(merge_vcf_path, div_vcf_path, output_vcf_path):
    div_vcf = VariantFile(div_vcf_path)
    merge_vcf = VariantFile(merge_vcf_path)
    # going to use the header from the div_vcf -- I'm not sure if in general one is preferable over the other
    # but in the benchmarks so far it seems they'll usually be the same
    hd = div_vcf.header
    out_vcf = VariantFile(output_vcf_path, 'w', header=hd)
    for vcf in [div_vcf, merge_vcf]:
        for rec in vcf.fetch():
            out_vcf.write(rec)
    out_vcf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Correct event positions in vcf recording div_dDUPs')
    parser.add_argument('--div_dDUP_vcf', help='Input vcf from R1 generation step')
    parser.add_argument('--merge_input_vcf', help='Input vcf with previously simulated events to merge into div. vcf')
    parser.add_argument('--merge_output_vcf', help='Output path for merged div. dDUP / other type vcf')
    parser.add_argument('--label_entire_event', action='store_true', help='Boolean flag for whether or not'
                                                                          'output vcf should report one record '
                                                                          'per div. repeat with start/end giving '
                                                                          'the overall span of the entire event ('
                                                                          'rather than interval A or a '
                                                                          'specifically)')
    args = parser.parse_args()

    incorporate_events(args.merge_input_vcf, args.div_dDUP_vcf, args.merge_output_vcf)
    correct_positions(args.merge_output_vcf, args.label_entire_event)
