from pysam import VariantFile
import argparse
import numpy as np
from bisect import bisect


def correct_positions_div(input_vcf, label_entire_event=False, avoid_intervals=False):
    # reads input vcf, adjusts POS, END, TARGET values to account for previous events in the reference
    in_vcf = VariantFile(input_vcf)
    hd = in_vcf.header
    file_suffix = '_avoid_intervals.vcf' if avoid_intervals else '_pos_corrected.vcf'
    out_vcf = VariantFile(input_vcf[:-4] + file_suffix, 'w', header=hd)
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
        # logic by which records are output dictated by the kind of output vcf we want
        if not avoid_intervals:
            if label_entire_event and rec.id == 'div_dDUP':
                # simplified record just reporting start and end positions with respect to entire "A_A*" region
                rec.stop = rec.info['TARGET'] + rec.info['SVLEN']
                rec.info['TARGET'] = -1
                rec.info['SVLEN'] = rec.stop - rec.start
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
                rec_A.info['DIV_REPEAT'] = 'NA'
                # setting the above unsets rec.stop, so need to reset it to the original
                rec_A.stop = rec.stop
                rec_B.id = 'div_dDUP_B'
                rec_B.alts = ('div_dDUP_B',)
                rec_B.start = rec_B.info['TARGET']
                rec_B.stop = rec_B.info['TARGET'] + rec.info['SVLEN']
                rec_B.info['SVTYPE'] = 'div_dDUP_B'
                rec_B.info['TARGET'] = -1
                rec_B.info['DIV_REPEAT'] = 'NA'
                out_vcf.write(rec_A)
                out_vcf.write(rec_B)

    out_vcf.close()
    in_vcf.close()


# TODO: I think this can be deleted, I'm pretty sure I was mistaken that you'll have to adjust positions in this way..
def correct_positions_orig_events(orig_vcf_path, r1_vcf_path):
    # when placing divergent repeats into an already edited reference, need to
    # adjust the R1 (reference with added repeats) positions to account for events
    # already added to the reference that will be edited and turned into R2
    orig_vcf = VariantFile(orig_vcf_path)
    r1_vcf = VariantFile(r1_vcf_path)
    hd = r1_vcf.header
    out_vcf = VariantFile(r1_vcf_path[:-4] + '_ORIGEVENTS_SHIFT.vcf', 'w', header=hd)
    orig_recs = {}
    for rec in orig_vcf.fetch():
        orig_recs[rec.start] = rec
    # walk through original recs and define a piecewise function that tells you how many places to
    # shift depending on where you are on the genome
    # ---> this will then be used to shift the R1 events by going to each event position and shifting by
    # f(position) for f the above shift function ^
    shift = 0
    r1_shifts = {}
    for i in range(len(list(orig_recs.keys()))):
        pos = sorted(orig_recs.keys())[i]
        r = orig_recs[pos]
        svlen = int(r.info['SVLEN'])
        if r.id == 'DEL':
            shift -= svlen
            r1_shifts[pos + svlen] = shift
        # TODO: THE POSITION ADJUSTMENT WONT WORK IF INSERTIONS AREN'T INCLUDED IN THE VCF (which I think they're not...)
        #  --> they're in the bed file that was output by the simulation -- rerun complex_bed_to_vcf.py including INS in the allow list
        elif r.id in ['DUP', 'INS']:
            shift += svlen
            r1_shifts[pos + svlen] = shift
        elif r.id == 'dDUP_A':
            # for dDUP_A, the shift will be of size len(dDUP_A), but will only be applied to events after the dDUP_B
            shift += svlen
            # need to look down the list and look for the adjacent dDUP_B (might lie beyond interceding simple events)
            for j in range(i, len(list(orig_recs.keys()))):
                # stop at the first dDUP_B we see
                if orig_recs[sorted(orig_recs.keys())[j]].id == 'dDUP_B':
                    r_ddup_b = orig_recs[sorted(orig_recs.keys())[j]]
                    # double check that we're conceivably next to the dDUP_A of concern
                    # *** this will raise an exception if a dDUP_A is unmatches
                    assert np.abs(r_ddup_b.start - r.stop) < 100, 'dDUP_A found without a mate'
                    r1_shifts[r_ddup_b.stop] = shift
                    break
    # walk through r1 records and apply a shift to each one determined by the r1_shifts dict
    # --> the right shift value will be the dict value matching the greatest key that is < the r1 record's position
    for rec in r1_vcf.fetch():
        supremum_index = bisect(sorted(r1_shifts.keys()), rec.start) - 1
        event_shift = r1_shifts[sorted(r1_shifts.keys())[supremum_index]]
        rec.start += event_shift
        # all recs in this edit step should be div_dDUPs
        if rec.id == 'div_dDUP':
            rec.info['TARGET'] += event_shift
        out_vcf.write(rec)



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
