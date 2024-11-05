#!/usr/bin/env python3

"""insilicoSV: simulator of structural variants.

(redesigned version)
"""

import argparse
from collections import defaultdict, Counter
import copy
import functools
import logging
import os
import os.path
import random
import shutil
import time
from typing_extensions import Any, Optional, cast

from pysam import FastaFile
import yaml

from insilicosv import utils
from insilicosv.utils import (
    Locus, Region, RegionSet, RegionFilter, OverlapMode, chk, error_context,
    n_valid_placements, has_duplicates, if_not_none)
from insilicosv.sv_defs import (SV, Operation, Transform, TransformType,
                                Breakend, BreakendRegion)
from insilicosv.variant_set_makers import make_variant_set_from_config
from insilicosv.output import OutputWriter

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class SVSimulator:

    """
    Simulates SVs from a given config.
    """

    config: dict[str, Any]
    svs: list[SV]
    rois: RegionSet
    rois_list: list[Region]
    rois_list_pos: int
    used_sv_regions: RegionSet
    translocated_chroms: set[str]
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

        self.convert_legacy_config()
        self.pre_check_config()

        self.verbose = self.config['sim_settings'].get('verbose', False)
        if self.verbose:
            loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
            for a_logger in loggers:
                a_logger.setLevel(logging.DEBUG)

        if 'random_seed' in self.config['sim_settings']:
            random_seed = self.config['sim_settings']['random_seed']
            chk(isinstance(random_seed, (type(None), int, float, str)),
                'invalid random_seed')
            if random_seed == 'time':
                random_seed = time.time_ns() % 1000000000
            logger.info(f'random seed: {random_seed}')
            random.seed(random_seed)


        self.output_path = os.path.dirname(config_path)
        
        self.reference = FastaFile(self.config['sim_settings']['reference'])
        self.chrom_lengths = { chrom: chrom_length
                               for chrom, chrom_length in zip(self.reference.references,
                                                              self.reference.lengths) }

        self.used_sv_regions = RegionSet()
        self.translocated_chroms = set()
        
    def pre_check_config(self) -> None:
        config = self.config

        chk(isinstance(config, dict), 'Config must be a dict')

        for k in config:
            chk(k in ('sim_settings', 'variant_sets', 'overlap_regions', 'blacklist_regions'),
                f'invalid top-level config: {k}')

        chk('sim_settings' in config, 'Missing sim_settings')
        sim_settings = config['sim_settings']
        chk(isinstance(sim_settings, dict), 'sim_settings must be a dict')

        for k in sim_settings:
            chk(k in ('reference', 'max_tries', 'homozygous_only', 'min_intersv_dist',
                      'random_seed', 'output_no_haps', 'output_paf', 'output_svops_bed',
                      'output_paf_intersv', 'verbose'),
                f'invalid sim_settings key: {k}')

        chk('reference' in sim_settings, 'Missing reference under sim_settings')
        chk(utils.is_readable_file(sim_settings['reference']), f'reference must be a readable file')
        ref_idx = sim_settings['reference'] + '.fai'
        chk(utils.is_readable_file(ref_idx), f'reference index must be a readable file')
        chk(isinstance(config.get('variant_sets'), list),
            'variant_sets must specify a list')

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
                    chk(not key.startswith('overlap_'),
                        f'Using {key} in {variant_set} requires specifying overlap_regions in global config')

    def convert_legacy_config(self) -> None:
        chk(isinstance(self.config, dict), 'Config must be a dict')
        chk('sim_settings' in self.config, 'Missing sim_settings')

        for k in ('prioritize_top', 'fail_if_placement_issues'):
            if k in self.config['sim_settings']:
                del self.config['sim_settings'][k]
                logger.warning(f'Please remove deprecated setting {k}')

        for variant_set in self.config.get('variant_sets', []):
            variant_set_orig = copy.deepcopy(variant_set)
            with error_context(variant_set):
                cfg = variant_set
                chk(isinstance(cfg, dict), 'each variant set must be a dict')
                if 'min_length' in cfg and 'max_length' in cfg:
                    chk('length_ranges' not in cfg, 'cannot combine old and new syntax for length ranges')
                    chk(utils.is_list_of((int, type(None)), cfg['min_length']), f'bad min_length')
                    chk(utils.is_list_of((int, type(None)), cfg['max_length']), f'bad max_length')
                    chk(len(cfg['min_length']) == len(cfg['max_length']),
                        f'min/max length lists must have same length')
                    cfg['length_ranges'] = [
                        [min_len, max_len] for min_len, max_len in zip(cfg['min_length'], cfg['max_length'])]
                    del cfg['min_length']
                    del cfg['max_length']
                    logger.warning(f'Please update to use length_ranges: {variant_set_orig}')

                for old, new in (('overlap_type', 'overlap_mode'),
                                 ('vcf_path', 'import')):
                    if old in cfg:
                        chk(new not in cfg, f'cannot combine old and new syntax for {new}')
                        cfg[new] = cfg[old]
                        del cfg[old]
                        logger.warning(f'Please update to use {new}: {variant_set_orig}')

                if cfg.get('type') == 'TRA_UNBALANCED':
                    cfg['type'] = 'TRA_NONRECIPROCAL'
                    logger.warning(f'Please update to use TRA_NONRECIPROCAL: {variant_set_orig}')
                if cfg.get('type') == 'TRA_BALANCED':
                    cfg['type'] = 'TRA_RECIPROCAL'
                    logger.warning(f'Please update to use TRA_RECIPROCAL: {variant_set_orig}')


    def load_rois(self) -> None:
        self.rois = RegionSet.from_beds(
            utils.as_list(self.config.get('overlap_regions', [])),
            verbose=self.verbose)

        # we model unconstrained placement as constrained placement with
        # overlap mode 'contained' and roi 'all reference'
        all_reference = (
            RegionSet.from_fasta(self.config['sim_settings']['reference'],
                                 region_kind='_reference_'))
        self.rois.add_region_set(all_reference)

        self.rois_list = list(self.rois.get_region_list())
        random.shuffle(self.rois_list)
        self.rois_list_pos = 0

        self.blacklist_regions = RegionSet()
        for blacklist_region_file in utils.as_list(self.config.get('blacklist_regions', [])):
            if blacklist_region_file.lower().endswith('.bed'):
                self.blacklist_regions.add_region_set(RegionSet.from_beds([blacklist_region_file]))
            elif blacklist_region_file.lower().endswith('.vcf'):
                self.blacklist_regions.add_region_set(RegionSet.from_vcf(blacklist_region_file))
            else:
                chk(f'Cannot import blacklist regions from {blacklist_region_file}: '
                    f'unsupported file type')

    def run(self) -> None:

        self.construct_svs()
        self.place_svs()
        self.output_results()

    def construct_svs(self) -> None:
        self.svs = []
        for vset_num, variant_set_config in enumerate(self.config['variant_sets']):
            vset_svs = make_variant_set_from_config(variant_set_config, self.config['sim_settings'])
            for sv in vset_svs:
                sv.info['VSET'] = vset_num
            self.svs.extend(vset_svs)

        assert not has_duplicates(sv.sv_id for sv in self.svs)

    def place_svs(self) -> None:

        self.load_rois()

        # TODO: globally optimize SV placement
        
        # Currently, we place SVs one at a time.

        random.shuffle(self.svs)

        self.determine_sv_placement_order()

        t_start_placing = time.time()
        logger.info(f'Placing {len(self.svs)} svs')
        t_last_status = time.time()
        for sv_num, sv in enumerate(self.svs):
            t_start_placing_sv = time.time()
            logger.debug(f'Placing {sv_num=} {sv=}')
            with error_context(sv.config_descr):
                self.place_sv(sv)
            logger.debug(f'Placed {sv_num=} {sv=} in {time.time()-t_start_placing_sv}s')
            assert sv.is_placed()
            for operation in sv.operations:
                if operation.chromosomal_translocation_source_breakend is not None:
                    assert operation.source_region is not None
                    self.translocated_chroms |= {operation.source_region.chrom,
                                                 operation.target_region.chrom}

            for region in sv.get_regions():
                region_padded = region.padded(self.config['sim_settings'].get('min_intersv_dist', 1))

                self.used_sv_regions.add_region(region_padded)
                new_roi_pieces = self.rois.chop(region_padded)

                for new_roi_piece in new_roi_pieces:
                    self.rois_list.insert(random.randint(0, len(self.rois_list)), new_roi_piece)

            if time.time() - t_last_status > 10:
                logger.info(f'Placed {sv_num} of {len(self.svs)} SVs in {time.time()-t_start_placing:.1f}s')
                t_last_status = time.time()
        logger.info(f'Placed {len(self.svs)} svs in {time.time()-t_start_placing:.1f}s.')

    def determine_sv_placement_order(self) -> None:
        # place most constrained SVs first
        logger.debug(f'Deciding placement order for {len(self.svs)} SVs')
        relevant_rois: dict[RegionFilter, list[Region]] = {}
        for sv in self.svs:
            with error_context(sv.config_descr):
                if sv.fixed_placement:
                    sv.num_valid_placements = 1
                else:
                    assert sv.roi_filter is not None
                    if sv.roi_filter not in relevant_rois:
                        logger.debug(f'getting relevant rois for {sv.roi_filter}')
                        relevant_rois[sv.roi_filter] = self.rois.filtered(
                            region_filter=sv.roi_filter).get_region_list()
                        logger.debug(f'got relevant rois for {sv.roi_filter}')

                    assert sv.overlap_mode is not None

                    sv.num_valid_placements = self.estimate_num_valid_placements(
                        regions=relevant_rois[sv.roi_filter],
                        overlap_mode=sv.overlap_mode,
                        anchor_length=sv.get_anchor_length())
                chk(sv.num_valid_placements > 0, f'No valid placements for SV')

        self.svs.sort(key=lambda sv: sv.num_valid_placements)

    def is_placement_valid(self, sv: SV, placement: list[Locus],
                           roi: Optional[Region] = None,
                           blacklist_regions: Optional[RegionSet] = None) -> bool:
        """Returns False if proposed placement of `sv` would run off chromosome or 
        touch blacklisted or already-used regions."""
        is_chromosomal_translocation = any(operation.chromosomal_translocation_source_breakend is not None
                                           for operation in sv.operations)
        for op_region in sv.get_regions(placement, roi):
            if self.used_sv_regions.overlaps_region(op_region):
                return False

            if blacklist_regions is not None and blacklist_regions.overlaps_region(op_region):
                return False

            if op_region.start < 0 or op_region.end > self.chrom_lengths[op_region.chrom]:
                return False

            if utils.percent_N(self.reference.fetch(reference=op_region.chrom,
                                                    start=op_region.start,
                                                    end=op_region.end)) > 0.05:
                return False

            if is_chromosomal_translocation and op_region.chrom in self.translocated_chroms:
                return False

        return True
    # end: def is_placement_valid(...)

    @functools.cache
    def get_relevant_blacklist_regions(self, blacklist_filter: Optional[RegionFilter]) -> RegionSet:
        return (RegionSet() if blacklist_filter is None else
                self.blacklist_regions.filtered(region_filter=blacklist_filter))

    @staticmethod
    def estimate_num_valid_placements(regions: list[Region], anchor_length: Optional[int],
                                      overlap_mode: OverlapMode) -> int:
        if overlap_mode == OverlapMode.EXACT:
            return len(regions)
        assert anchor_length is not None
        n_plc = 0
        for region in regions:
            n_plc += n_valid_placements(region_length=region.length(),
                                        anchor_length=anchor_length,
                                        overlap_mode=overlap_mode)
            if n_plc > 10000:
                break
        return n_plc

    def choose_anchor_placement(
            self, roi_filter: Optional[RegionFilter],
            anchor_length: Optional[int],
            overlap_mode: OverlapMode) -> tuple[Optional[Region], Optional[Region]]:
        """Finds the next ROI that meets `roi_filter` (if given) and on which at least
        one anchor placement of `anchor_length` satisfying `overlap_mode` is possible,
        randomly chooses an anchor placement on the ROI, and returns the pair (anchor, roi).
        If no suitable placement exists for any ROI, returns (None, None).
        """

        n_rois_tried = 0

        while n_rois_tried < len(self.rois_list):
            self.rois_list_pos = (self.rois_list_pos + 1) % len(self.rois_list)
            n_rois_tried += 1

            roi = self.rois_list[self.rois_list_pos]
            if roi not in self.rois:
                continue

            if roi_filter and not roi_filter.satisfied_for(roi):
                continue

            assert (0 <= roi.orig_start <= roi.start <= roi.end <= roi.orig_end
                    <= self.chrom_lengths[roi.chrom])

            if overlap_mode == OverlapMode.EXACT:
                assert anchor_length is None
                if (roi.start != roi.orig_start or
                    roi.end != roi.orig_end):
                    continue
                return (roi, roi)

            assert anchor_length is not None

            n_plc = n_valid_placements(anchor_length=anchor_length, overlap_mode=overlap_mode,
                                       region_length=roi.length())
            if not n_plc:
                continue

            if overlap_mode == OverlapMode.CONTAINED:
                offset = random.randint(0, n_plc-1)
                if anchor_length == 0:
                    offset += 1
                anchor_region = roi.replace(start=roi.start + offset,
                                            end=roi.start + offset + anchor_length)
                return (anchor_region, roi)

            if overlap_mode == OverlapMode.PARTIAL:
                assert anchor_length > 1
                assert roi.length() > 1
                max_overlap = min(anchor_length-1, roi.length()-1)
                assert max_overlap > 0

                if (roi.start != roi.orig_start and 
                    roi.end != roi.orig_end):
                    continue

                offset = random.randint(1, max_overlap)
                anchor_region_choices = []
                if roi.start == roi.orig_start:
                    anchor_region_choices.append(roi.replace(start=roi.start + offset - anchor_length,
                                                             end=roi.start + offset))
                if roi.end == roi.orig_end:
                    anchor_region_choices.append(roi.replace(start=roi.end - offset,
                                                             end=roi.end - offset + anchor_length))
                anchor_region = random.choice(anchor_region_choices)
                return (anchor_region, roi)
        # end: while self.rois_list_pos != start_rois_list_pos:

        return (None, None)

    # end: def choose_anchor_placement(...)

    def place_sv(self, sv: SV) -> None:

        assert not sv.is_placed()

        if sv.fixed_placement:
            chk(self.is_placement_valid(sv, sv.fixed_placement),
                f'cannot place imported SV')
            sv.set_placement(placement=sv.fixed_placement, roi=None)
            return

        assert sv.overlap_mode is not None and sv.anchor is not None

        n_placement_attempts = 0
        max_tries = self.config["sim_settings"].get("max_tries", 1000)

        blacklist_regions = self.get_relevant_blacklist_regions(sv.blacklist_filter)

        while not sv.is_placed() and (n_placement_attempts < max_tries):
            n_placement_attempts += 1
            
            anchor_region, anchor_roi = (
                self.choose_anchor_placement(
                    roi_filter=sv.roi_filter,
                    anchor_length=sv.get_anchor_length(),
                    overlap_mode=sv.overlap_mode))

            chk(anchor_region is not None, f"No anchor placements for variant")
            assert anchor_region is not None and anchor_roi is not None
            assert anchor_region.chrom == anchor_roi.chrom

            breakend_interval_lengths = list(sv.breakend_interval_lengths)
            if sv.overlap_mode == OverlapMode.EXACT:
                assert (breakend_interval_lengths[sv.anchor.start_breakend] is None and
                        sv.breakend_interval_min_lengths[sv.anchor.start_breakend] is None)
                breakend_interval_lengths[sv.anchor.start_breakend] = anchor_region.length()
                
            placement_dict: dict[Breakend, Locus] = {
                sv.anchor.start_breakend: Locus(chrom=anchor_region.chrom, pos=anchor_region.start)}

            def place_interchrom(breakend_locus: Locus) -> Locus:
                new_chrom_choices = [chrom for chrom, chrom_length in self.chrom_lengths.items()
                                     if chrom != breakend_locus.chrom and chrom_length >= 2]
                chk(new_chrom_choices, f'interchromosomal events require >1 chrom of length >1')
                new_chrom = random.choice(new_chrom_choices)
                return Locus(chrom=new_chrom,
                             pos=random.randint(1, self.chrom_lengths[new_chrom]-1))

            breakend_locus = placement_dict[sv.anchor.start_breakend]
            for breakend in range(sv.anchor.start_breakend-1, -1, -1):

                if breakend_interval_lengths[breakend] is None and sv.is_interchromosomal:
                    breakend_locus = place_interchrom(breakend_locus)
                else:
                    if breakend_interval_lengths[breakend] is not None:
                        breakend_interval_length = breakend_interval_lengths[breakend]
                    else:
                        # unbounded intra-chromosomal dispersion
                        min_dispersion_len = if_not_none(sv.breakend_interval_min_lengths[breakend], 0)
                        breakend_interval_length = random.randint(
                            min_dispersion_len, max(min_dispersion_len, breakend_locus.pos))
                    assert isinstance(breakend_interval_length, int)
                    breakend_locus = breakend_locus.shifted(-breakend_interval_length)
                placement_dict[breakend] = breakend_locus

            breakend_locus = placement_dict[sv.anchor.start_breakend]
            for breakend in range(sv.anchor.start_breakend+1, len(breakend_interval_lengths)+1):
                if breakend_interval_lengths[breakend-1] is None and sv.is_interchromosomal:
                    breakend_locus = place_interchrom(breakend_locus)
                else:
                    if breakend_interval_lengths[breakend-1] is not None:
                        breakend_interval_length = breakend_interval_lengths[breakend-1]
                    else:
                        min_dispersion_len = if_not_none(sv.breakend_interval_min_lengths[breakend-1], 0)
                        breakend_interval_length = random.randint(
                            min_dispersion_len,
                            max(min_dispersion_len,
                                self.chrom_lengths[breakend_locus.chrom] - breakend_locus.pos))
                    assert isinstance(breakend_interval_length, int)
                    breakend_locus = breakend_locus.shifted(breakend_interval_length)
                placement_dict[breakend] = breakend_locus

            placement: list[Locus] = [placement_dict[breakend]
                                      for breakend in sv.breakends]

            if self.is_placement_valid(sv, placement, anchor_roi, blacklist_regions):
                sv.set_placement(placement=placement, roi=anchor_roi)
        # end: while not sv.is_placed() and (n_placement_attempts < max_tries):

        if not sv.is_placed():
            raise RuntimeError(f'Could not place SV within {max_tries=}')

        if sv.is_placed():
            logger.debug(f'Placed {sv=} within {n_placement_attempts=}')

    # end: def place_sv(self, sv)

    def output_results(self) -> None:
        logger.info('Writing outputs')
        output_writer = OutputWriter(self.svs, self.reference, self.chrom_lengths,
                                     self.output_path, self.config)
        output_writer.output_vcf()
        output_writer.output_haps()
        output_writer.output_novel_insertions()
        output_writer.output_stats()

        # dev only
        output_writer.output_svops_bed()

        self.reference.close()

        logger.info(f'Output path: {self.output_path}')

    # for testing only
    def produce_variant_genome(self, fasta1_out, fasta2_out, ins_fasta, bedfile, stats_file=None, initial_reset=True,
                               verbose=False, export_to_file=True):
        self.run()
        shutil.copyfile(os.path.join(self.output_path, 'sim.hapA.fa'), fasta1_out)
        shutil.copyfile(os.path.join(self.output_path, 'sim.hapB.fa'), fasta2_out)
        shutil.copyfile(os.path.join(self.output_path, 'sim.novel_insertions.fa'), ins_fasta)

# end: class SVSimulator

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--config', required=True, help='yaml config file')
    
    return parser.parse_args()


def run_simulator():
    start_time = time.time()

    logger.info('insilicoSV version 1.0')

    args = parse_args()

    simulator = SVSimulator(config_path=args.config)
    
    simulator.run()

    logger.info(f'insilicoSV finished in {time.time()-start_time:.1f}s')

if __name__ == '__main__':
    run_simulator()
