from abc import ABC, abstractmethod
from dataclasses import dataclass
from copy import copy
from enum import Enum
from functools import cached_property
import logging
from typing_extensions import TypeAlias, Optional, Any, cast

from insilicosv.utils import (
    Locus, Region, RegionSet, OverlapMode, RegionFilter, chk, pairwise)

class TransformType(Enum):

    IDENTITY = "IDENTITY"
    INV = "INV"
    DEL = "DEL"

@dataclass(frozen=True)
class Transform:
    transform_type: TransformType
    is_in_place: bool

    n_copies: int = 1
    divergence_prob: float = 0
    replacement_seq: Optional[str] = None

    def __post_init__(self) -> None:
        assert (self.transform_type != TransformType.DEL or
                (self.is_in_place and self.divergence_prob == 0 and self.replacement_seq is None))
        assert 0 <= self.divergence_prob <= 1
        assert not (self.divergence_prob > 0 and self.replacement_seq is not None)
        assert self.n_copies >= 0
        assert not self.is_in_place or self.n_copies == 1

Breakend: TypeAlias = int

@dataclass(frozen=True)
class BreakendRegion:

    # includes both its ends

    start_breakend: Breakend
    end_breakend: Breakend

    def __post_init__(self) -> None:
        assert 0 <= self.start_breakend <= self.end_breakend

@dataclass
class Operation:
    transform: Transform

    source_breakend_region: Optional[BreakendRegion] = None
    novel_insertion_seq: Optional[str] = None
    chromosomal_translocation_source_breakend: Optional[Breakend] = None
    chromosomal_translocation_active: bool = False

    target_insertion_breakend: Optional[Breakend] = None
    target_insertion_order: Optional[tuple] = None

    placement: Optional[list[Locus]] = None

    op_info: Optional[dict] = None

    @property
    def transform_type(self):
        return self.transform.transform_type

    @property
    def is_in_place(self):
        return self.transform.is_in_place

    def assert_valid(self):
        assert ((self.source_breakend_region is not None) + (self.novel_insertion_seq is not None) +
                (self.chromosomal_translocation_source_breakend is not None)) == 1
        assert (self.novel_insertion_seq is None and self.chromosomal_translocation_source_breakend is None) or not self.is_in_place
        assert (not self.is_in_place) == (self.target_insertion_breakend is not None and
                                          self.target_insertion_order is not None)
        assert (self.source_breakend_region is None or
                self.source_breakend_region.end_breakend == self.source_breakend_region.start_breakend + 1)

    def get_source_region(self, placement: list[Locus]) -> Optional[Region]:
        if self.source_breakend_region is not None:
            assert (placement[self.source_breakend_region.start_breakend].chrom ==
                    placement[self.source_breakend_region.end_breakend].chrom)
            return Region(chrom=placement[self.source_breakend_region.start_breakend].chrom,
                          start=placement[self.source_breakend_region.start_breakend].pos,
                          end=placement[self.source_breakend_region.end_breakend].pos)
        if self.chromosomal_translocation_source_breakend is not None:
            return Region(chrom=placement[self.chromosomal_translocation_source_breakend].chrom,
                          start=placement[self.chromosomal_translocation_source_breakend].pos,
                          end=placement[self.chromosomal_translocation_source_breakend].pos)
        return None

    def get_target_region(self, placement: list[Locus]) -> Region:
        if self.is_in_place:
            return cast(Region, self.get_source_region(placement))
        else:
            assert self.target_insertion_breakend is not None
            assert self.target_insertion_order is not None
            return Region(chrom=placement[self.target_insertion_breakend].chrom,
                          start=placement[self.target_insertion_breakend].pos,
                          end=placement[self.target_insertion_breakend].pos,
                          order_key=self.target_insertion_order)

    @cached_property
    def source_region(self) -> Optional[Region]:
        assert self.placement is not None
        return self.get_source_region(self.placement)

    @cached_property
    def target_region(self) -> Region:
        assert self.placement is not None
        return self.get_target_region(self.placement)

# end: class Operation

@dataclass
class SV(ABC):
    # an ID for this SV, unique within this simulation
    sv_id: str

    # breakend_interval_lengths[i] is the fixed distance between breakends i and i+1,
    # or None if the distance is not fixed (for unbounded dispersions or for
    # components constrained to exactly coincide with an ROI).
    # the number of breakends of this SV is len(breakend_interval_lengths)+1.
    breakend_interval_lengths: list[Optional[int]]
    
    # minimum breakend interval lengths (for unbounded dispersions with lower length bounds)
    breakend_interval_min_lengths: list[Optional[int]]

    # whether unbounded dispersions in this SV are interchromosomal
    is_interchromosomal: bool

    # list of operations comprising this SV, specified in terms of its breakends.
    operations: list[Operation]

    #
    # Placement constraints
    #

    # criteria for selecting a ROI on which to place this SV
    roi_filter: Optional[RegionFilter]

    # sub-region of this SV which is constrained to overlap with the ROI
    anchor: Optional[BreakendRegion]

    # how `anchor` should overlap with the ROI
    overlap_mode: Optional[OverlapMode]

    # criteria for selecting blacklist regions with which no part of this SV may overlap
    blacklist_filter: Optional[RegionFilter]

    # a specific placement prescribed for this SV -- used for SVs imported from a vcf
    fixed_placement: Optional[list[Locus]]

    #########################

    # truth data about this SV, for the VCF INFO field
    info: dict[str, Any]

    # for each of two haplotypes, whether this SV is present on that haplotype
    genotype: tuple[bool, bool]

    # user-readable string identifying the variant set configuration from which
    # this SV was created (for error messages only)
    config_descr: str

    #
    # Fields set when SV is placed
    #

    placement: Optional[list[Locus]] = None
    roi: Optional[Region] = None

    # Fields used while finding a placement
    num_valid_placements: int = 0

    def __post_init__(self) -> None:
        assert self.sv_id
        assert len(self.breakend_interval_min_lengths) == len(self.breakend_interval_lengths)
        assert all(bi_min_len is None or bi_len is None
                   for bi_len, bi_min_len in zip(self.breakend_interval_lengths,
                                                 self.breakend_interval_min_lengths))
        for operation in self.operations:
            operation.assert_valid()
            assert (not operation.source_breakend_region or 
                    0 <= operation.source_breakend_region.start_breakend < len(self.breakend_interval_lengths))
        assert ((self.overlap_mode is None) == (self.anchor is None) == 
                (self.roi_filter is None))
        if self.anchor is not None:
            self.get_anchor_length()
            assert (
                0 <= self.anchor.start_breakend <= self.anchor.end_breakend 
                <= len(self.breakend_interval_lengths))
            assert (self.overlap_mode != OverlapMode.EXACT or
                    self.anchor.end_breakend == self.anchor.start_breakend + 1)
            assert (self.overlap_mode != OverlapMode.PARTIAL or
                    self.anchor.end_breakend > self.anchor.start_breakend)
        assert self.fixed_placement is None or self.overlap_mode is None
        assert self.genotype and sum(self.genotype)
        
        for operation in self.operations:
            if operation.target_insertion_order is not None:
                operation.target_insertion_order = (self.sv_id,) + operation.target_insertion_order

    @property
    def breakends(self) -> tuple[Breakend, ...]:
        return tuple(range(len(self.breakend_interval_lengths) + 1))

    @abstractmethod
    def to_vcf_records(self, sim_settings) -> list[dict]:
        raise NotImplementedError()

    def get_anchor_length(self) -> Optional[int]:

        assert self.anchor is not None
        if self.overlap_mode == OverlapMode.EXACT:
            chk(self.anchor.end_breakend == self.anchor.start_breakend + 1,
                f'overlap_mode "exact" requires that anchor constrain '
                f'exactly one symbol: {self.config_descr}')
            chk((self.breakend_interval_lengths[self.anchor.start_breakend] is None and
                 self.breakend_interval_min_lengths[self.anchor.start_breakend] is None),
                f'overlap_mode "exact" requires leaving the length '
                f'of the anchor symbol unspecified: {self.config_descr}')
            return None
        
        lengths_within_anchor = self.breakend_interval_lengths[
            self.anchor.start_breakend:self.anchor.end_breakend]
        chk(all(length is not None for length in lengths_within_anchor),
            f'All symbols within anchor must have a specified length: {self.config_descr}')
        anchor_length = sum(cast(list[int], lengths_within_anchor))
        chk(self.overlap_mode != OverlapMode.PARTIAL or anchor_length > 0,
            f'overlap_mode "partial" requires non-empty anchor: {self.config_descr}')

        return anchor_length

    def set_placement(self, placement: list[Locus], roi: Optional[Region]) -> None:
        placement = copy(placement)
        assert len(placement) == len(self.breakend_interval_lengths) + 1
        for breakend1, breakend2 in pairwise(self.breakends):
            assert (placement[breakend1].chrom != placement[breakend2].chrom or
                    placement[breakend1].pos <= placement[breakend2].pos)
            assert (self.breakend_interval_lengths[breakend1] is None or
                    (placement[breakend2].chrom == placement[breakend1].chrom and
                     placement[breakend2].pos - placement[breakend1].pos
                     == self.breakend_interval_lengths[breakend1]))
            assert (self.breakend_interval_min_lengths[breakend1] is None or
                    (placement[breakend2].chrom == placement[breakend1].chrom and
                     placement[breakend2].pos - placement[breakend1].pos
                     >= self.breakend_interval_min_lengths[breakend1]))

        self.placement = placement
        self.roi = roi

        for operation in self.operations:
            assert operation.placement is None
            operation.placement = placement

    def is_placed(self):
        return self.placement is not None

    def get_regions(self, placement: Optional[list[Locus]] = None, roi: Optional[Region] = None) -> list[Region]:
        if placement is None:
            placement = self.placement
        assert placement is not None
        if roi is None:
            roi = self.roi

        regions = []

        for operation in self.operations:
            source_region = operation.get_source_region(placement)
            target_region = operation.get_target_region(placement)
            if source_region is not None:
                regions.append(source_region)
            if (target_region is not None and
                target_region != source_region):
                regions.append(target_region)

        if self.overlap_mode == OverlapMode.EXACT:
            assert roi is not None
            regions.append(roi)

        regions_merged: list[Region] = []
        for region in sorted(set(regions)):
            if (regions_merged and
                region.chrom == regions_merged[-1].chrom and
                region.start == regions_merged[-1].end):
                regions_merged[-1] = regions_merged[-1].replace(end=region.end)
            else:
                regions_merged.append(region)

        return regions_merged

# end: class SV

