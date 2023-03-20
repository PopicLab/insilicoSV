from pysam import VariantFile
import argparse
from collections import defaultdict


def correct_positions_div(input_vcf, label_entire_event=False, avoid_intervals=False):
    # reads input vcf, adjusts POS, END, TARGET values to account for previous events in the reference
    in_vcf = VariantFile(input_vcf)
    hd = in_vcf.header
    file_suffix = '_decoy_divrepeat_intervals.vcf' if avoid_intervals else '_pos_corrected.vcf'
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
                # *** shifting the start position of a vcf record automatically shifts the stop position as well
                # --> when simulating further events on top of a div. repeat simulation, the current logic is such
                # --> that the second round of events will have their locations determined already accounting for the
                # --> insertions into the R2 reference -- that is, we won't need to shift their event positions
                evt.start += total_ins_len
                evt.info['TARGET'] += total_ins_len
                total_ins_len += evt.info['SVLEN']

    for chrom in vcf_recs.keys():
        for rec in vcf_recs[chrom].values():
            # logic by which records are output dictated by the kind of output vcf we want
            if not avoid_intervals:
                if label_entire_event and rec.id == 'div_dDUP':
                    # simplified record just reporting start and end positions with respect to entire "A_A*" region
                    rec.stop = rec.info['TARGET'] + rec.info['SVLEN']
                    rec.info['TARGET'] = -1
                    rec.info['SVLEN'] = rec.stop - rec.start
                # # remove divergent repeat sequence for reporting
                # rec.info['DIV_REPEAT'] = '(removed)'
                out_vcf.write(rec)
            else:
                # output separate records for each A and A* ("B") to serve as avoid intervals for further simulation
                if rec.id == 'div_dDUP':
                    rec_A = rec.copy()
                    rec_B = rec.copy()
                    rec_A.id = 'div_dDUP_A'
                    rec_A.alts = ('div_dDUP_A',)
                    rec_A.info['SVTYPE'] = 'div_dDUP_A'
                    rec_A.info['TARGET'] = -1
                    # rec_A.info['DIV_REPEAT'] = 'NA'
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
                    # rec_B.info['DIV_REPEAT'] = 'NA'
                    out_vcf.write(rec_A)
                    out_vcf.write(rec_B)

    out_vcf.close()
    in_vcf.close()


def incorporate_events(merge_vcf_path, div_vcf_path, output_vcf_path):
    div_vcf = VariantFile(div_vcf_path)
    merge_vcf = VariantFile(merge_vcf_path)
    hd = div_vcf.header
    out_vcf = VariantFile(output_vcf_path, 'w', header=hd)
    for vcf in [div_vcf, merge_vcf]:
        for rec in vcf.fetch():
            out_vcf.write(rec)
    out_vcf.close()
    div_vcf.close()
    merge_vcf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Correct event positions in vcf recording div_dDUPs')
    parser.add_argument('--div_dDUP_vcf', help='Input vcf from R1 generation step')
    parser.add_argument('--merge_input_vcf', default=None, help='Input vcf with previously simulated events to merge into div. vcf')
    parser.add_argument('--merge_output_vcf', default=None, help='Output path for merged div. dDUP / other type vcf')
    parser.add_argument('--label_entire_event', action='store_true', help='Boolean flag for whether or not'
                                                                          'output vcf should report one record '
                                                                          'per div. repeat with start/end giving '
                                                                          'the overall span of the entire event ('
                                                                          'rather than interval A or a '
                                                                          'specifically)')
    parser.add_argument('--avoid_intervals', action='store_true', help='Boolean flag to trigger the preparation of a'
                                                                       'vcf giving div_dDUP interval A and A* in'
                                                                       'separate records (for use as avoid_intervals'
                                                                       'vcf when simulating additional events after)')
    args = parser.parse_args()
    assert not (args.label_entire_event and args.avoid_intervals), 'Cannot provide label_entire_event and avoid_intervals' \
                                                                   'flags both as True'

    if args.merge_input_vcf is not None and args.merge_output_vcf is not None:
        incorporate_events(args.merge_input_vcf, args.div_dDUP_vcf, args.merge_output_vcf)
        correct_positions_div(args.merge_output_vcf, args.label_entire_event, args.avoid_intervals)
    else:
        correct_positions_div(args.div_dDUP_vcf, args.label_entire_event, args.avoid_intervals)
