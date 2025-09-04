#!/usr/bin/env python3

"""insilicoSV: simulator of structural variants.
"""

import argparse
import functools
import logging
import os
import os.path
import random
import shutil
import time
from typing_extensions import Any
import numpy as np
from intervaltree import Interval, IntervalTree
from pysam import FastaFile
import yaml
import copy
from collections import defaultdict
from math import floor

from insilicosv import utils, __version__
from insilicosv.utils import (
    Locus, Region, RegionSet, OverlapMode, chk,
    has_duplicates, if_not_none, pairwise)
from insilicosv.sv_defs import SV, Breakend, Transform, TransformType, Operation, BreakendRegion
from insilicosv.variant_set import make_variant_set_from_config
from insilicosv.output import OutputWriter

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

DEFAULT_PERCENT_N = 0.05
DEFAULT_MAX_TRIES = 100
MIN_INTERSV_DIST = 1
FILTER_SMALL_CHR = 0


class SVSimulator:
    """
    Simulates SVs from a given config.
    """

    config: dict[str, Any]
    svs: list[SV]
    # Dictionary which keys are variant set idx that have overlap constrained and containing a list of ROIs
    rois_overlap: dict[list[Region]]
    # Region set defined by the reference updated to keep track of the available regions.
    reference_regions: RegionSet
    # Region Set defined by the reference regions to keep track of available regions for overlap SVs
    reference_sv_overlap_regions: RegionSet
    blacklist_regions: RegionSet
    reference: FastaFile
    chrom_lengths: dict[str, int]
    output_path: str
    verbose: bool

    def __init__(self, config_path):
        try:
            with open(config_path) as config_yaml:
                self.config = yaml.safe_load(config_yaml)
        except Exception as exception:
            chk(False, f'Error reading config file {config_path}: {exception}')

        self.pre_check_config()

        self.verbose = self.config.get('verbose', False)
        if self.verbose:
            loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
            for a_logger in loggers:
                a_logger.setLevel(logging.DEBUG)
        if 'random_seed' in self.config:
            random_seed = self.config['random_seed']
            chk(isinstance(random_seed, (type(None), int, float, str)), 'invalid random_seed')
            if random_seed == 'time':
                random_seed = time.time_ns() % 1000000000
            logger.info(f'random seed: {random_seed}')
            random.seed(random_seed)
            np.random.seed(random_seed)
        # Variables used to prefilter ROIs according to the constraints of each category of SV
        self.overlap_ranges = {}
        self.overlap_kinds = {}
        self.overlap_modes = {}
        self.num_svs = {}
        self.rois_overlap = {}
        self.output_path = os.path.dirname(config_path)
        self.allow_hap_overlap = self.config.get('allow_hap_overlap', False)

        # Tree to keep track of the overlap SV regions with the operation performed to give to the OutputWriter
        self.has_overlap_sv = False
        self.overlap_sv_regions = RegionSet()

        self.reference = FastaFile(self.config['reference'])
        self.chrom_lengths = {chrom: chrom_length
                              for chrom, chrom_length in zip(self.reference.references,
                                                             self.reference.lengths)}

    def pre_check_config(self):
        config = self.config

        chk(isinstance(config, dict), 'Config must be a dict')

        for k in config:
            chk(k in ('reference', 'max_tries', 'max_random_breakend_tries', 'homozygous_only', 'heterozygous_only',
                      'min_intersv_dist',
                      'random_seed', 'output_no_haps', 'output_adjacencies', 'output_paf', 'output_svops_bed',
                      'th_proportion_N',
                      'output_paf_intersv', 'verbose', 'variant_sets', 'overlap_regions', 'blacklist_regions',
                      'filter_small_chr', 'allow_hap_overlap'),
                f'invalid top-level config: {k}', error_type='syntax')

        chk(utils.is_readable_file(config['reference']), f'reference must be a readable file')
        chk(isinstance(config.get('variant_sets'), list), 'variant_sets must specify a list')

        if 'overlap_regions' in config:
            config['overlap_regions'] = utils.as_list(config['overlap_regions'])
            chk(isinstance(config['overlap_regions'], list),
                'overlap_regions must be a list')
            chk(all(utils.is_readable_file(overlap_spec) for overlap_spec in config['overlap_regions']),
                'Each item in overlap_regions must be a readable file')
        else:
            for variant_set in config['variant_sets']:
                chk(isinstance(variant_set, dict), f'variant set must be a dict: {variant_set}')
                for key in variant_set.keys():
                    chk(not key.startswith('overlap_') or key == 'overlap_sv' or
                        variant_set.get('overlap_mode', None) in ['terminal', 'whole-chromosome'],
                        f'Using {key} in {variant_set} requires specifying overlap_regions in global config')

    def load_rois(self):
        self.reference_regions = RegionSet.from_fasta(self.config['reference'],
                                                      self.config.get('filter_small_chr', FILTER_SMALL_CHR),
                                                      region_kind='_reference_',
                                                      allow_hap_overlap=self.allow_hap_overlap)
        if self.has_overlap_sv:
            self.reference_sv_overlap_regions = copy.deepcopy(self.reference_regions)
        # Get the ROIs for overlap constraints
        if any(mode is not None for mode in self.overlap_modes.values()):
            rois_overlap = RegionSet.from_beds(utils.as_list(self.config.get('overlap_regions', [])),
                                               to_region_set=False)
            min_bounds = [min_bound for min_bound, _ in self.overlap_ranges.values()]
            global_min_bound = 0
            logger.info(
                f'Filtering out ROIs not satisfying the overlap constraints from initial {len(rois_overlap)} ROIs')
            if None not in min_bounds:
                global_min_bound = min(min_bounds)
                rois_overlap = sorted(rois_overlap, key=lambda x: x.length(), reverse=True)
            # Get the reference regions
            n_removed_rois = 0
            for roi_index, roi in enumerate(rois_overlap):
                # Check if the roi is under the global minimum, in which case no following ROI would satisfy the constraints.
                if roi.length() < global_min_bound:
                    n_removed_rois += len(rois_overlap) - roi_index
                    break
                if not roi.chrom in self.reference_regions.chrom2itree:
                    n_removed_rois += 1
                    continue
                added_roi = False
                for sv_idx, sv_category in enumerate(self.overlap_ranges):
                    if self.overlap_modes[sv_idx] in [None, OverlapMode.TERMINAL, OverlapMode.CHROM]: continue
                    if roi.length() < if_not_none(self.overlap_ranges[sv_category][0], 0): continue
                    if (self.overlap_modes[sv_category] in [OverlapMode.CONTAINING, OverlapMode.EXACT] and
                            roi.length() > if_not_none(self.overlap_ranges[sv_category][1], roi.length() + 1)):
                        continue
                    if not roi.kind in self.overlap_kinds[sv_category] and not 'all' in self.overlap_kinds[sv_category]:
                        found = False
                        for kind in self.overlap_kinds[sv_category]:
                            if kind in roi.kind:
                                found = True
                                break
                        if not found: continue
                    added_roi = True
                    self.rois_overlap[sv_category].append(roi)
                if not added_roi: n_removed_rois += 1

            for sv_idx, sv_category in enumerate(self.overlap_ranges):
                if self.overlap_modes[sv_idx] not in [OverlapMode.TERMINAL, OverlapMode.CHROM]: continue
                for chrom, chrom_length in self.chrom_lengths.items():
                    self.rois_overlap[sv_category].append(Region(chrom=chrom, start=0, end=chrom_length, kind='chr',
                                                                 orig_start=0, orig_end=chrom_length))

            for sv_category in self.rois_overlap:
                # Check if there is enough ROIs to fit all the SVs of one category independently
                error_message_num_rois = ("Only {} ROIs satisfying the constraints "
                                          "(overlap mode: {}, type of ROIs: {}, overlap range: {}) of the variant_set {} containing {} SVs").format(
                    len(self.rois_overlap[sv_category]), self.overlap_modes[sv_category],
                    self.overlap_kinds[sv_category],
                    self.overlap_ranges[sv_category], sv_category, self.num_svs[sv_category])

                hap_overlap_mult = 2 if self.allow_hap_overlap else 1
                if self.overlap_modes[sv_category] in [OverlapMode.CONTAINING, OverlapMode.EXACT]:
                    # We have one containing and exact overlap per ROI and haplotype
                    chk(self.num_svs[sv_category] <= hap_overlap_mult * len(self.rois_overlap[sv_category]),
                        error_message_num_rois)
                elif self.overlap_modes[sv_category] == OverlapMode.PARTIAL:
                    # We can have up tow partial SVs per ROI and per haplotype.
                    chk(self.num_svs[sv_category] <= hap_overlap_mult * 2 * len(self.rois_overlap[sv_category]),
                        error_message_num_rois)
                elif self.overlap_modes[sv_category] == OverlapMode.CHROM:
                    # We have one CHROM overlap per chromosome and haplotype
                    chk(self.num_svs[sv_category] <= hap_overlap_mult * len(self.rois_overlap[sv_category]),
                        error_message_num_rois)
                elif self.overlap_modes[sv_category] == OverlapMode.TERMINAL:
                    # We have two TERMINAL overlaps per chromosome and haplotype
                    chk(self.num_svs[sv_category] <= 2 * hap_overlap_mult * len(self.rois_overlap[sv_category]),
                        error_message_num_rois)
                # Shuffle the ROIs so the selection is not biased on their positions in the input bed file
                random.shuffle(self.rois_overlap[sv_category])
            logger.info(f'{n_removed_rois} ROIs filtered')
        self.blacklist_regions = RegionSet()
        for blacklist_region_file in utils.as_list(self.config.get('blacklist_regions', [])):
            logger.info(f'Processing blacklist region file {blacklist_region_file}')
            if blacklist_region_file.lower().endswith('.bed'):
                self.blacklist_regions.add_region_set(RegionSet.from_beds([blacklist_region_file], to_region_set=True))
            elif blacklist_region_file.lower().endswith('.vcf'):
                self.blacklist_regions.add_region_set(RegionSet.from_vcf(blacklist_region_file))
            else:
                chk(f'Cannot import blacklist regions from {blacklist_region_file}: '
                    f'unsupported file type, please provide a .bed or .vcf file', error_type='type')
            logger.info(f'Blacklist region file {blacklist_region_file} processed.')

    def run(self):
        self.construct_svs()
        self.load_rois()
        self.place_svs()
        self.output_results()

    def construct_svs(self):
        self.svs = []
        logger.info('Constructing SVs from {} categories'.format(len(self.config['variant_sets'])))
        for vset_num, variant_set_config in enumerate(self.config['variant_sets']):
            variant_set_config['VSET'] = vset_num
            vset_svs, ranges, kinds, mode, header = make_variant_set_from_config(variant_set_config, self.config)
            for sv in vset_svs:
                sv.info['VSET'] = vset_num
            self.overlap_ranges[vset_num] = ranges
            self.overlap_kinds[vset_num] = kinds
            self.overlap_modes[vset_num] = mode
            self.svs.extend(vset_svs)

            if vset_svs[0].allow_sv_overlap:
                self.has_overlap_sv = True

            self.num_svs[vset_num] = len(vset_svs)
        logger.info(f'Constructed {len(self.svs)} SVs')
        self.rois_overlap = {vset_num: [] for vset_num in range(len(self.config['variant_sets']))}
        assert not has_duplicates(sv.sv_id for sv in self.svs)

    def update_available_reference(self, sv):
        for region in sv.get_regions():
            region_padded = region.padded(self.config.get('min_intersv_dist', MIN_INTERSV_DIST))
            self.reference_regions.chop(region_padded, sv.genotype)

            if sv.overlap_mode == OverlapMode.CHROM and sv.info['OP_TYPE'] == 'DUP':
                orig_op = sv.operations[0]
                operations = []
                for copy_num in range(orig_op.transform.n_copies):
                    # For the writing of the output, we create an operation per copy to create new chromosome copies
                    transform = Transform(
                        transform_type=TransformType.IDENTITY,
                        is_in_place=False,
                        divergence_prob=orig_op.transform.divergence_prob,
                        n_copies=1
                    )
                    operations.append(
                        Operation(transform=transform,
                                  op_info=orig_op.op_info,
                                  source_breakend_region=BreakendRegion(start_breakend=Breakend(0),
                                                                        end_breakend=Breakend(1)),
                                  target_insertion_breakend=Breakend(2),
                                  placement=[Locus(chrom=region.chrom, pos=region.start),
                                             Locus(chrom=region.chrom, pos=region.end),
                                             Locus(
                                                 chrom=region.chrom + f'_copy_{copy_num}',
                                                 pos=region.start)]
                                  )
                    )
                sv.operations = operations

    def update_overlap_svs(self, sv):
        for region in sv.get_regions():
            self.overlap_sv_regions.add_region(region, sv=sv, allow_hap_overlap=self.allow_hap_overlap)
            self.reference_sv_overlap_regions.chop(region, sv.genotype)

    def place_svs(self):
        # Currently, we place SVs one at a time.
        self.determine_sv_placement_order()
        t_start_placing = time.time()
        logger.info(f'Placing {len(self.svs)} svs')
        t_last_status = time.time()
        # Starting ROI index for the categories with containing or exact overlaps.
        roi_indexes = [0 for _ in self.rois_overlap]
        for sv_num, sv in enumerate(self.svs):
            t_start_placing_sv = time.time()
            logger.debug(f'Placing {sv_num=} {sv=}')
            svset = sv.info['VSET']
            roi_indexes[svset] = self.place_sv(sv, roi_indexes[svset])
            logger.debug(f'Placed {sv_num=} {sv=} in {time.time() - t_start_placing_sv}s')
            assert sv.is_placed()

            if not sv.allow_sv_overlap:
                self.update_available_reference(sv)
            else:
                self.update_overlap_svs(sv)

            if time.time() - t_last_status > 10:
                logger.info(f'Placed {sv_num} of {len(self.svs)} SVs in {time.time() - t_start_placing:.1f}s')
                t_last_status = time.time()
        logger.info(f'Placed {len(self.svs)} svs in {time.time() - t_start_placing:.1f}s.')

    def determine_sv_placement_order(self) -> None:
        # place most constrained SVs first
        logger.info(f'Deciding placement order for {len(self.svs)} SVs')
        types_order = ['SV OVERLAP', 'FIXED', OverlapMode.CHROM, OverlapMode.TERMINAL, OverlapMode.EXACT,
                       OverlapMode.PARTIAL, OverlapMode.CONTAINING, OverlapMode.CONTAINED, None]
        for sv in self.svs:
            if sv.allow_sv_overlap:
                sv.priority = 0
            elif sv.fixed_placement:
                sv.priority = 1
            else:
                distance = sum([dist for dist in sv.breakend_interval_lengths if dist is not None]) + 2
                sv.priority = (types_order.index(sv.overlap_mode)) + 1 / distance
        self.svs.sort(key=lambda sv: sv.priority)

    def is_placement_valid(self, sv, placement):
        """Returns False if proposed placement of `sv` would run off chromosome or
        touch blacklisted regions."""
        # The placement does not set all breakend positions
        if not len(placement) == len(sv.breakend_interval_lengths) + 1: return False
        chk(all([locus.pos <= self.chrom_lengths[locus.chrom] for locus in placement]),
            'Please make sure that the imported'
            ' SV positions are within the chromosome length,'
            f' provided {sv}', error_type='value')
        for breakend1, breakend2 in pairwise(sv.breakends):
            # Ensure the positions are ordered in a same chromosome
            if not (placement[breakend1].chrom != placement[breakend2].chrom or
                    placement[breakend1].pos <= placement[breakend2].pos): return False
            # Ensure the distances between breakends are satisfied
            if not (sv.breakend_interval_lengths[breakend1] is None or
                    (placement[breakend2].chrom == placement[breakend1].chrom and
                     placement[breakend2].pos - placement[breakend1].pos
                     == sv.breakend_interval_lengths[breakend1])): return False
            # Ensure the minimum distance between breakends is satisfied
            if not (sv.breakend_interval_min_lengths[breakend1] is None or
                    (placement[breakend2].chrom == placement[breakend1].chrom and
                     placement[breakend2].pos - placement[breakend1].pos
                     >= sv.breakend_interval_min_lengths[breakend1])): return False
        for op_region in sv.get_regions(placement):
            # Ensure the regions covered by the SV do not contain a proportion of Ns above th_proportion_N
            # To ensure that am insertion target is not in between two Ns, the region is padded
            min_bound = max(0, op_region.start - 1)
            max_bound = min(self.chrom_lengths[op_region.chrom], op_region.end + 1)
            if utils.percent_N(self.reference.fetch(reference=op_region.chrom,
                                                    start=min_bound,
                                                    end=max_bound)) > self.config.get('th_proportion_N',
                                                                                      DEFAULT_PERCENT_N):
                return False
        return True

    # end: def is_placement_valid(...)

    @functools.cache
    def get_relevant_blacklist_regions(self, blacklist_filter):
        return (RegionSet() if blacklist_filter is None else
                self.blacklist_regions.filtered(region_filter=blacklist_filter))

    def get_breakend(self, hap_id, reference_regions, containing_region=None, avoid_chrom=None, blacklist_regions=None,
                     roi_length=0, total_length=0):
        max_random_tries = self.config.get("max_random_breakend_tries", DEFAULT_MAX_TRIES)
        num_tries = 0
        breakend = None
        ref_roi = None
        chromosomes = [chrom for chrom in self.reference_regions.chrom2itree if chrom != avoid_chrom]
        while breakend is None or ref_roi is None and num_tries < max_random_tries:
            breakend, ref_roi = self.get_random_breakend(reference_regions, containing_region=containing_region,
                                                         chromosomes=chromosomes,
                                                         blacklist_regions=blacklist_regions, roi_length=roi_length,
                                                         total_length=total_length, hap_id=hap_id)
            num_tries += 1

        if breakend is None:
            # The breakend failed to be assigned by random sampling, look for a region fitting the SV actually available.
            return self.get_breakend_from_regions(reference_regions, containing_region=containing_region,
                                                  avoid_chrom=avoid_chrom,
                                                  blacklist_regions=blacklist_regions, roi_length=roi_length,
                                                  total_length=total_length,
                                                  hap_id=hap_id)
        return breakend, ref_roi

    def get_random_breakend(self, reference_regions, containing_region=None, chromosomes=None, blacklist_regions=None,
                            roi_length=0, total_length=0, hap_id=0):
        if containing_region is None:
            chromosomes = [chrom for chrom in chromosomes if self.chrom_lengths[chrom] - total_length >= 0]
            bounds = {chrom: (0, self.chrom_lengths[chrom] - total_length) for chrom in chromosomes}
            chrom = chromosomes[random.randint(0, len(chromosomes) - 1)]
        else:
            chrom = containing_region.chrom
            bounds = {chrom: (containing_region.start, containing_region.end)}

        breakend = random.randint(bounds[chrom][0], bounds[chrom][1])
        region = Region(start=breakend, end=breakend, chrom=chrom)
        if (chrom in blacklist_regions.chrom2itree) and (
                blacklist_regions.strictly_contains_point(region.start, region.chrom) or
                blacklist_regions.strictly_contains_point(region.end, region.chrom)):
            return None, None

        ref_rois = self.get_reference_interval(region, reference_regions, hap_id)
        ref_roi = None
        for roi in ref_rois:
            # There might be two overlapping ref intervals if there was an insertion target and the min_intersv_dist is 0
            if breakend + roi_length > roi.data.end: continue
            ref_roi = roi.data
            break
        return region, ref_roi

    def get_breakend_from_regions(self, reference_regions, containing_region=None, avoid_chrom=None,
                                  blacklist_regions=None,
                                  roi_length=0, total_length=0, hap_id=0):
        max_random_tries = self.config.get("max_random_breakend_tries", DEFAULT_MAX_TRIES)
        # The region is not constrained, we use the interval tree defined from the reference file
        chrom_trees = reference_regions.chrom2itree
        rois = []
        weights = []
        for chrom_tree, tree in chrom_trees.items():
            if (avoid_chrom is not None) and (chrom_tree == avoid_chrom): continue
            ref_intervals = tree
            if containing_region is not None:
                # The breakend is after a dispersion and has to satisfy the constraints induced by already placed breakends.
                if chrom_tree != containing_region.chrom: continue
                ref_intervals = []
                overlaps = tree[hap_id].overlap(containing_region.start - 0.2, containing_region.end + 0.2)
                for overlap in overlaps:
                    # Chop the intervals overlapping the containing region
                    min_interval = min(overlap.data.start, containing_region.start)
                    max_interval = max(overlap.data.end, containing_region.end)
                    ref_intervals.append(Interval(begin=min_interval - 0.1, end=max_interval + 0.1,
                                                  data=overlap.data.replace(start=min_interval, end=max_interval)))
            for interval in ref_intervals:
                # Remove the intervals that are too small or too close to the end of the chromosome
                if interval.length() < roi_length: continue
                if interval.data.start + total_length > self.chrom_lengths[chrom_tree]: continue
                rois.append((interval, interval))
                weights.append(interval.length() + 1)
        if not rois:
            return None, None

        total_weights = sum(weights)
        invalid = True
        num_iteration = 0
        roi = None
        ref_roi = None
        breakend = None
        # Ensure the breakends of the anchor are not in a blacklist region
        while invalid and num_iteration < max_random_tries:
            roi_idx = np.random.choice(len(rois), p=[weight / total_weights for weight in weights])
            roi, ref_roi = rois[roi_idx]
            # We are looking for a breakend.
            max_bound = min(roi.data.end - roi_length, self.chrom_lengths[roi.data.chrom] - total_length)
            if roi.data.start > max_bound: continue

            if roi.data.start == max_bound:
                breakend = max_bound
            else:
                breakend = np.random.randint(roi.data.start, max_bound)

            invalid = False
            if blacklist_regions and roi.data.chrom in blacklist_regions.chrom2itree:
                invalid = blacklist_regions.strictly_contains_point(breakend, roi.chrom)
            num_iteration += 1
        if invalid:
            return None, None

        roi = Interval(begin=breakend - 0.1, end=breakend + 0.1, data=roi.data.replace(start=breakend, end=breakend))
        return roi.data, ref_roi.data

    # Provides a list of valid rois and weights corresponding to their length to uniformly draw from
    def get_overlap_region(self, sv_category, roi_index, init_roi, reference_regions, hap_id, anchor_length=None,
                           overlap_mode=None, roi_filter=None):
        # The region is constrained we use the interval tree defined from the bed file
        roi_list = self.rois_overlap[sv_category]
        # Add the beginning of the ROIs list at the end of the list of ROIs to check to ensure all are checked in case we run out.
        for region in (roi_list[roi_index:] + roi_list[:init_roi]):
            roi_index = (roi_index + 1) % len(roi_list)
            if not roi_filter.satisfied_for(region): continue
            valid_region, ref_roi = self.check_interval_overlap(region, reference_regions, roi_filter, anchor_length,
                                                                overlap_mode, hap_id)
            if valid_region is not None:
                # A copy of the valid region is returned so we do not change the provided overlap regions (several SVs might overlap it)
                return copy.deepcopy(valid_region), ref_roi.data, roi_index
        return None, None, None

    def check_interval_overlap(self, region, reference_regions, roi_filter, anchor_length, overlap_mode, hap_id):
        # Only keep regions of the interval that are in the reference
        ref_intervals = self.get_reference_interval(region, reference_regions, hap_id)

        random.shuffle(ref_intervals)
        if not ref_intervals:
            return None, None
        if overlap_mode in [OverlapMode.CONTAINED, OverlapMode.TERMINAL]:
            # Remove too small regions if the overlap is contained or terminal
            for ref_interval in ref_intervals:
                if region.length() <= anchor_length: continue
                left_bound = max(region.start, ref_interval.data.start)
                right_bound = min(region.end, ref_interval.data.end)
                intersection = region.replace(start=left_bound, end=right_bound)

                if intersection.length() < anchor_length: continue
                # Discards intervals smaller than the minimum overlap.
                if (roi_filter.region_length_range[0] is not None) and (
                        intersection.length() < roi_filter.region_length_range[0]): continue
                return intersection, ref_interval
            return None, None
        # If the overlap is exact or chrom we disregard the intervals that have already been sliced or are not fully covered by the reference
        elif overlap_mode in [OverlapMode.EXACT, OverlapMode.CHROM]:
            if ((ref_intervals is None) or (len(ref_intervals) != 1) or
                    (region.orig_start != region.start) or (region.orig_end != region.end)):
                return None, None

            ref_interval = ref_intervals.pop()
            if (ref_interval is None) or (ref_interval.end < region.end) or (ref_interval.begin > region.start):
                return None, None
            return region, ref_interval
        # If the overlap is partial we remove regions that have been sliced on both sides (only allows contained overlap)
        elif overlap_mode == OverlapMode.PARTIAL:
            # If the length of the anchor is 1 a strict overlap is not possible by definition
            if anchor_length == 1:
                return None, None
            # In case the regions defined by the bed file do not match the  ones defined by the reference we might have no or several containing intervals
            for ref_interval in ref_intervals:
                # Check that the region is within the limits of the reference
                left_bound = max(region.start, ref_interval.data.start)
                right_bound = min(region.end, ref_interval.data.end)
                if right_bound <= left_bound: continue
                # Check if the constraint region has an extremity overlapped by the reference, if not there is no possible partial overlap
                if (right_bound == ref_interval.data.end) and (left_bound == ref_interval.data.start): continue
                intersection = region.replace(start=left_bound, end=right_bound)
                # There is at least one valid partial overlap on the left or the right
                overlap_left = ref_interval.data.start + anchor_length < intersection.end - 1
                overlap_right = ref_interval.data.end - anchor_length > intersection.start + 1
                if not (overlap_left or overlap_right): continue
                # Discards regions smaller than the minimum overlap.
                if (roi_filter.region_length_range[0] is not None) and (
                        intersection.length() < roi_filter.region_length_range[0]): continue
                return intersection, ref_interval
            return None, None
        # If the overlap is containing remove the intervals that have been sliced or too big to be strictly included on both sides
        elif overlap_mode == OverlapMode.CONTAINING:
            if ((region.orig_start != region.start) or
                    (region.orig_end != region.end) or
                    region.length() > anchor_length - 2):
                return None, None
            # check that the interval is fully contained in the reference
            if len(ref_intervals) > 1: return None, None
            ref_interval = ref_intervals.pop()

            if (ref_interval.data.start > region.start) or (ref_interval.data.end < region.end):
                return None, None
            return region, ref_interval
        else:
            chk(False, f'Invalid overlap_mode {overlap_mode}')

    def get_reference_interval(self, roi, reference_regions, hap_id):
        # Get the reference interval overlapping a ROI if any.
        for chrom_tree, trees in reference_regions.chrom2itree.items():
            if chrom_tree == roi.chrom:
                # Padding to ensure that intervals reduced to a point are still processed correctly.
                overlap_interval = list(trees[hap_id].overlap(roi.start - 0.2, roi.end + 0.2))
                # For determinism as the overlap function output is a set
                overlap_interval.sort(key=lambda x: (x.begin, x.end))
                # To prevent biases by always selecting the first interval.
                random.shuffle(overlap_interval)
                return overlap_interval

    # From input roi and ref_roi, places the anchor in roi such that the overlap constraints are fulfilled  and
    # the anchor fits in ref_roi which represents a region non used by another SV.
    def choose_anchor_placement(self, roi, ref_roi, anchor_length, overlap_mode, region_length_range=[None, None],
                                blacklist_regions=None):
        """Finds the next ROI that meets `roi_filter` (if given) and on which at least
        one anchor placement of `anchor_length` satisfying `overlap_mode` is possible,
        randomly chooses an anchor placement on the ROI, and returns the pair (anchor, roi).
        If no suitable placement exists for any ROI, returns (None, None).
        """
        max_random_tries = self.config.get("max_random_breakend_tries", DEFAULT_MAX_TRIES)

        if overlap_mode in [OverlapMode.EXACT, OverlapMode.CHROM, OverlapMode.TERMINAL]:
            if overlap_mode == OverlapMode.TERMINAL:
                # Randomly choose an extremity of the arm
                chrom_length = self.chrom_lengths[roi.chrom]
                start = 0
                end = anchor_length
                if (random.randint(0, 1) and chrom_length == ref_roi.end) or start < ref_roi.start:
                    start = chrom_length - anchor_length
                    end = chrom_length
                roi = Region(chrom=roi.chrom, start=start, end=end)

            # Ensure the breakends of the roi are not in a blacklist region
            if blacklist_regions and roi.chrom in blacklist_regions.chrom2itree:
                invalid = blacklist_regions.strictly_contains_point(roi.start, roi.chrom)
                invalid += blacklist_regions.strictly_contains_point(roi.end, roi.chrom)
                if invalid:
                    return None, None
            return roi, ref_roi

        if overlap_mode == OverlapMode.CONTAINED:
            invalid = True
            num_iteration = 0
            anchor_region = None
            # Ensure the breakends of the anchor are not in a blacklist region
            while invalid and num_iteration < max_random_tries:
                bound_placement = roi.end - roi.start - anchor_length
                offset = random.randint(0, bound_placement)
                anchor_region = roi.replace(start=roi.start + offset,
                                            end=roi.start + offset + anchor_length)
                invalid = False
                if blacklist_regions and anchor_region.chrom in blacklist_regions.chrom2itree:
                    invalid = blacklist_regions.strictly_contains_point(anchor_region.start, anchor_region.chrom)
                    invalid += blacklist_regions.strictly_contains_point(anchor_region.end, anchor_region.chrom)
                num_iteration += 1
            if invalid:
                return None, None
            return anchor_region, ref_roi

        if overlap_mode == OverlapMode.PARTIAL:
            # Looking for the possible positions for the start of the anchor
            possible_starts = []
            # Overlap on the left side
            # The minimum requested overlap has to be satisfied if provided.
            left_overlap = max(ref_roi.start, roi.start - anchor_length + if_not_none(region_length_range[0], 1))
            # On the right we ensure that the anchor won't end after the end of the region in which case it would be a containing overlap.
            # The maximum overlap length requested has to be satisfied if provided.
            right_overlap = min(roi.start - if_not_none(region_length_range[1], 1), roi.end - anchor_length - 1)
            if (roi.start == roi.orig_start) and (right_overlap >= left_overlap):
                possible_starts += list(range(left_overlap, right_overlap + 1))

            # Overlap on the right side.
            left_overlap = max(roi.end - if_not_none(region_length_range[1], anchor_length) + 1,
                               roi.start + 1)
            right_overlap = min(ref_roi.end - anchor_length, roi.end - if_not_none(region_length_range[0], 1))
            if (roi.end == roi.orig_end) and (right_overlap >= left_overlap):
                possible_starts += list(range(left_overlap, right_overlap + 1))

            # Prevent exact overlap
            if roi.start in possible_starts and anchor_length == roi.length():
                possible_starts.remove(roi.start)

            if not possible_starts: return (None, None)

            invalid = True
            num_iteration = 0
            anchor_region = None
            # Ensure the breakends of the anchor are not in a blacklist region
            while invalid and num_iteration < max_random_tries:
                start_breakend = random.choice(possible_starts)
                end_breakend = start_breakend + anchor_length

                anchor_region = roi.replace(start=start_breakend, end=end_breakend)

                invalid = False
                if blacklist_regions and anchor_region.chrom in blacklist_regions.chrom2itree:
                    invalid = blacklist_regions.strictly_contains_point(anchor_region.start, anchor_region.chrom)
                    invalid += blacklist_regions.strictly_contains_point(anchor_region.end, anchor_region.chrom)
                num_iteration += 1
            if invalid:
                return None, None
            return anchor_region, ref_roi

        if overlap_mode == OverlapMode.CONTAINING:
            left_bound = max(ref_roi.start, roi.end + 1 - anchor_length)
            right_bound = min(roi.start - 1, ref_roi.end - anchor_length)

            if left_bound > right_bound: return None, None

            invalid = True
            num_iteration = 0
            anchor_region = None
            # Ensure the breakends of the anchor are not in a blacklist region
            while invalid and num_iteration < max_random_tries:
                start_breakend = random.randint(left_bound, right_bound)
                end_breakend = start_breakend + anchor_length

                if end_breakend > ref_roi.end: return None, None

                anchor_region = roi.replace(start=start_breakend, end=end_breakend)

                invalid = False
                if blacklist_regions and anchor_region.chrom in blacklist_regions.chrom2itree:
                    invalid = blacklist_regions.strictly_contains_point(anchor_region.start, anchor_region.chrom)
                    invalid += blacklist_regions.strictly_contains_point(anchor_region.end, anchor_region.chrom)
                num_iteration += 1
            if invalid:
                return None, None
            return anchor_region, ref_roi

    # Sum the SV lengths up to the next dispersion
    def sum_lengths(self, breakend_interval_lengths, breakends, dispersions):
        contiguous_length = 0
        for idx, breakend in enumerate(breakends):
            if (breakend in dispersions) or (breakend_interval_lengths[idx] is None): break
            contiguous_length += breakend_interval_lengths[idx]
        return contiguous_length

    # Function to traverse the breakends starting from an anchor region or the first breakend. Breakends are placed using a
    # known distance or a new breakend is randomly chosen among the available regions. When starting from an anchor,
    # the traversal can be in forward or backward order.
    def propagate_placement(self, placement_dict, roi, ref_roi, breakends, lengths, min_dist, blacklist_regions,
                            interchromosomal, dispersions, backward, hap_id, anchor_breakends=None):
        shift = 1 if not backward else -1
        # In case of exact overlap with several symbols in the anchor
        anchor_roi = roi
        for pos, breakend in enumerate(breakends):
            if breakend + shift in placement_dict:
                # If we have an anchor so we do not move its breakend.
                locus = placement_dict[breakend + shift]
                roi = Region(chrom=locus.chrom, start=locus.pos, end=locus.pos, kind=roi.kind, motif=roi.motif)
                continue
            distance = lengths[pos]
            if distance is None:
                # If we have an interchromosomal dispersion we want to change chromosome otherwise keep the same one.
                avoid_chrom = roi.chrom
                containing_region = None
                in_anchor = (anchor_breakends is not None and
                             anchor_breakends.start_breakend <= breakend < anchor_breakends.end_breakend)
                contiguous_length = self.sum_lengths(lengths[pos:], breakends[breakend:], dispersions[breakend:])
                # If interchromosomal the length is the contiguous length
                total_length = contiguous_length
                if not interchromosomal or in_anchor:
                    total_length = sum([length for length in lengths[pos:] if length is not None])
                    avoid_chrom = None
                    bound = if_not_none(min_dist[pos], 0)
                    if in_anchor:
                        # we add a breakend position between the last breakend placed from th same SV, and the end of the anchor roi
                        # +1 to ensure the previous symbol in the anchor doesn't have a length of 0
                        left_bound = roi.start + 1
                        # -1 to ensure the current symbol doesn't have a length of 0
                        right_bound = anchor_roi.end - 1
                    else:
                        left_bound = 0 if backward else (roi.end + bound)
                        right_bound = (roi.start - bound) if backward else self.chrom_lengths[roi.chrom]
                    if left_bound > right_bound: return None
                    containing_region = Region(chrom=roi.chrom, start=left_bound, end=right_bound, kind=roi.kind,
                                               motif=roi.motif)

                roi, ref_roi = self.get_breakend(reference_regions=self.reference_regions,
                                                 avoid_chrom=avoid_chrom,
                                                 blacklist_regions=blacklist_regions,
                                                 containing_region=containing_region,
                                                 roi_length=contiguous_length,
                                                 total_length=total_length,
                                                 hap_id=hap_id)
                if roi is None: return None
                position = roi.start
            else:
                position = roi.start - distance if backward else roi.start + distance
                if (blacklist_regions is not None) and (roi.chrom in blacklist_regions.chrom2itree):
                    # blacklist regions are the same on both haplotypes
                    if blacklist_regions.strictly_contains_point(position, roi.chrom):
                        return None
                if breakend in dispersions:
                    # The placement has to be valid on the haplotypes corresponding to the SV's genotype
                    overlap = self.reference_regions.chrom2itree[roi.chrom][hap_id].overlap(
                        Interval(begin=position - 0.2, end=position + 0.2))
                    valid = len(overlap) > 0
                    if valid:
                        ref_roi = overlap.pop().data
                else:
                    valid = ref_roi.start <= position <= ref_roi.end
                if not valid: return None
                roi = Region(start=position, end=position, chrom=roi.chrom)

            placement_dict[breakend + shift] = Locus(chrom=roi.chrom, pos=position)
        return placement_dict

    def place_sv(self, sv, roi_index):
        assert not sv.is_placed()
        reference_regions = self.reference_regions
        if sv.allow_sv_overlap:
            reference_regions = self.reference_sv_overlap_regions

        if sv.fixed_placement:
            chk(self.is_placement_valid(sv, sv.fixed_placement),
                f'cannot place imported SV {sv}, please check your SVs '
                f'are non overlapping and try lowering the min_intersv_dist or'
                f' increasing th_proportion_N.')
            sv.set_placement(placement=sv.fixed_placement, roi=None)
            return

        n_placement_attempts = 0
        max_tries = self.config.get("max_tries", DEFAULT_MAX_TRIES)
        blacklist_regions = self.get_relevant_blacklist_regions(sv.blacklist_filter)
        sv_set = sv.info['VSET']
        hap_id = 0
        if self.allow_hap_overlap:
            hap_id = 2 if sv.genotype[0] and sv.genotype[1] else sv.genotype[1]
        init_roi = 0
        if sv.overlap_mode in [OverlapMode.CONTAINED, OverlapMode.PARTIAL]:
            # Where to start checking the ROIs, prevent the bias of checking the first ROIs over and over
            roi_index = init_roi = random.randint(0, len(self.rois_overlap[sv_set]))

        while not sv.is_placed() and (n_placement_attempts < max_tries):
            n_placement_attempts += 1
            placement_dict: dict[Breakend, Locus] = {}
            breakend_interval_lengths = list(sv.breakend_interval_lengths)
            min_distances = list(sv.breakend_interval_min_lengths)
            anchor_start = anchor_end = 0
            if sv.anchor is not None:
                roi, ref_roi, roi_index = self.get_overlap_region(sv_category=sv_set,
                                                                  anchor_length=sv.get_anchor_length(),
                                                                  reference_regions=reference_regions,
                                                                  overlap_mode=sv.overlap_mode,
                                                                  roi_filter=sv.roi_filter,
                                                                  hap_id=hap_id,
                                                                  roi_index=roi_index,
                                                                  init_roi=init_roi)

                if roi is None or ref_roi is None:
                    chk(False, f'No available ROI satisfying the constraints for {sv}' +
                        f' of anchor length {sv.get_anchor_length()}' * (sv.get_anchor_length() is not None))
                anchor_start = sv.anchor.start_breakend
                anchor_end = sv.anchor.end_breakend
                roi, ref_roi = (
                    self.choose_anchor_placement(
                        roi=roi,
                        ref_roi=ref_roi,
                        anchor_length=sv.get_anchor_length(),
                        overlap_mode=sv.overlap_mode,
                        region_length_range=sv.roi_filter.region_length_range,
                        blacklist_regions=blacklist_regions))

                if (roi is None) or (ref_roi is None):
                    continue
            else:
                # Compute the length needed in the ROI to fit the breakends not seperated by dispersions
                total_length = contiguous_length = self.sum_lengths(breakend_interval_lengths,
                                                                    range(len(breakend_interval_lengths)),
                                                                    sv.dispersions)
                if not sv.is_interchromosomal:
                    total_length = sum([length for length in breakend_interval_lengths if length is not None])
                roi, ref_roi = self.get_breakend(reference_regions=reference_regions, hap_id=hap_id,
                                                 blacklist_regions=blacklist_regions,
                                                 roi_length=contiguous_length, total_length=total_length)
                if roi is None or ref_roi is None: break

            placement_dict[anchor_start] = Locus(chrom=roi.chrom, pos=roi.start)
            placement_dict[anchor_end] = Locus(chrom=roi.chrom, pos=roi.end)
            # Place the other breakends using the distances starting from the anchor ones towards the extremities.
            if anchor_start > 0:
                placement_dict = self.propagate_placement(hap_id=hap_id,
                                                          placement_dict=placement_dict,
                                                          roi=roi,
                                                          ref_roi=ref_roi,
                                                          breakends=[i for i in range(anchor_start, 0, -1)],
                                                          lengths=[breakend_interval_lengths[i] for i in
                                                                   range(anchor_start - 1, -1, -1)],
                                                          min_dist=[min_distances[i] for i in
                                                                    range(anchor_start - 1, -1, -1)],
                                                          blacklist_regions=blacklist_regions,
                                                          interchromosomal=sv.is_interchromosomal,
                                                          dispersions=sv.dispersions,
                                                          backward=True)
                if placement_dict is None: continue
            placement_dict = self.propagate_placement(hap_id=hap_id,
                                                      placement_dict=placement_dict,
                                                      roi=roi,
                                                      ref_roi=ref_roi,
                                                      breakends=[i for i in
                                                                 range(anchor_start, len(breakend_interval_lengths))],
                                                      lengths=breakend_interval_lengths[anchor_start:],
                                                      min_dist=min_distances[anchor_start:],
                                                      blacklist_regions=blacklist_regions,
                                                      interchromosomal=sv.is_interchromosomal,
                                                      dispersions=sv.dispersions,
                                                      backward=False,
                                                      anchor_breakends=sv.anchor)
            if placement_dict is None: continue
            placement: list[Locus] = [placement_dict[breakend] for breakend in sv.breakends]
            if not self.is_placement_valid(sv, placement): continue
            sv.set_placement(placement=placement, roi=roi, operation=sv.operations[0])
        # end: while not sv.is_placed() and (n_placement_attempts < max_tries):
        if not sv.is_placed():
            raise RuntimeError(f'Could not place SV within {max_tries=}')

        if sv.is_placed():
            logger.debug(f'Placed {sv=} within {n_placement_attempts=}')
        return roi_index

    # end: def place_sv(self, sv)

    def output_results(self) -> None:
        logger.info('Writing outputs')
        output_writer = OutputWriter(self.svs, self.overlap_sv_regions, self.reference, self.chrom_lengths,
                                     self.output_path, self.allow_hap_overlap, self.config)
        logger.info('Writing new haplotypes')
        output_writer.output_haps()
        logger.info('Writing VCF file')
        output_writer.output_vcf()
        logger.info('Writing novel insertions file')
        output_writer.output_novel_insertions()
        logger.info('Writing novel adjacencies file')
        output_writer.output_novel_adjacencies()
        output_writer.output_stats()

        # dev only
        output_writer.output_svops_bed()

        self.reference.close()

        logger.info(f'Output path: {self.output_path}')

    # for testing only
    def produce_variant_genome(self, fasta1_out, fasta2_out, ins_fasta):
        self.run()
        shutil.copyfile(os.path.join(self.output_path, 'sim.hapA.fa'), fasta1_out)
        shutil.copyfile(os.path.join(self.output_path, 'sim.hapB.fa'), fasta2_out)
        shutil.copyfile(os.path.join(self.output_path, 'sim.novel_insertions.fa'), ins_fasta)


# end: class SVSimulator

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--config', required=True, help='YAML config file')

    return parser.parse_args()


def run_simulator():
    start_time = time.time()
    logger.info('insilicoSV version %s' % __version__)
    args = parse_args()
    simulator = SVSimulator(config_path=args.config)
    simulator.run()
    logger.info(f'insilicoSV finished in {time.time() - start_time:.1f}s')


if __name__ == '__main__':
    run_simulator()