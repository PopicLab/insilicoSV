from abc import ABC, abstractmethod
from dataclasses import dataclass
from copy import copy
from enum import Enum
from functools import cached_property
from typing_extensions import TypeAlias, Optional, Any, cast, override

from insilicosv.utils import (
    Locus, Region, OverlapMode, RegionFilter, chk)

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

    def __post_init__(self):
        assert (self.transform_type != TransformType.DEL or
                (self.is_in_place and self.divergence_prob == 0 and self.replacement_seq is None))
        assert 0 <= self.divergence_prob <= 1
        assert self.divergence_prob == 0 or self.replacement_seq is None

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

    def get_source_region(self, placement):
        if self.source_breakend_region is not None:
            return Region(chrom=placement[self.source_breakend_region.start_breakend].chrom,
                          start=placement[self.source_breakend_region.start_breakend].pos,
                          end=placement[self.source_breakend_region.end_breakend].pos)
        return None

    def get_target_region(self, placement):
        if self.is_in_place:
            return cast(Region, self.get_source_region(placement))
        else:
            return Region(chrom=placement[self.target_insertion_breakend].chrom,
                          start=placement[self.target_insertion_breakend].pos,
                          end=placement[self.target_insertion_breakend].pos,
                          order_key=self.target_insertion_order)

    @cached_property
    def source_region(self):
        assert self.placement is not None
        return self.get_source_region(self.placement)

    @cached_property
    def target_region(self):
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

    # Positions of the dispersions
    dispersions: list[Optional[int]]

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

    # how `anchors` should overlap with the ROI
    overlap_mode: Optional[OverlapMode]

    # criteria for selecting blacklist regions with which no part of this SV may overlap
    blacklist_filter: Optional[RegionFilter]

    # a specific placement prescribed for this SV -- used for SVs imported from a vcf
    fixed_placement: Optional[Locus]

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

    def __post_init__(self):
        assert self.sv_id
        assert len(self.breakend_interval_min_lengths) == len(self.breakend_interval_lengths)
        assert all(bi_min_len is None or bi_len is None
                   for bi_len, bi_min_len in zip(self.breakend_interval_lengths,
                                                 self.breakend_interval_min_lengths))

        if self.overlap_mode == OverlapMode.CONTAINED:
            chk((self.roi_filter.region_length_range[0] is None) or (self.anchor.length() > self.roi_filter.region_length_range[0]),
                f'The anchor length is smaller than the minimum overlap for a contained overlap.')
            chk((self.roi_filter.region_length_range[1] is None) or (self.anchor.length() < self.roi_filter.region_length_range[1]),
                f'The anchor length is larger than the maximum overlap for a contained overlap.')
        if self.overlap_mode == OverlapMode.CONTAINING:
            chk((self.roi_filter.region_length_range[0] is None) or (self.anchor.length() > self.roi_filter.region_length_range[0]),
                f'The anchor length is smaller than the minimum overlap for a containing overlap.')
        if self.overlap_mode == OverlapMode.CONTAINING:
            chk((self.roi_filter.region_length_range[0] is None) or (self.anchor.length() >= self.roi_filter.region_length_range[0]),
                f'The anchor length is smaller than the minimum overlap for a partial overlap.')

        if self.anchor is not None:
            chk(all(
                length is not None for idx, length in enumerate(self.breakend_interval_lengths) if
                (idx not in self.dispersions) and
                ((self.overlap_mode != OverlapMode.EXACT) or not (
                            self.anchor.start_breakend <= idx < self.anchor.end_breakend))),
                f'A length range can onl be [null, null] for dispersions or the anchor for an exact overlap.')
            chk(all([breakend not in self.dispersions for breakend in range(self.anchor.start_breakend, self.anchor.end_breakend-1)]),
                f'anchors cannot contain a dispersion: {self.config_descr}')
            chk((self.overlap_mode != OverlapMode.EXACT) or
                all((self.breakend_interval_lengths[breakend] is None and
                 self.breakend_interval_min_lengths[breakend] is None) for breakend in range(self.anchor.start_breakend, self.anchor.end_breakend)),
                f'overlap_mode "exact" requires leaving the length of the anchor symbols unspecified: {self}')
            chk((self.overlap_mode != OverlapMode.PARTIAL and self.overlap_mode != OverlapMode.CONTAINING)
                or self.anchor.end_breakend != self.anchor.start_breakend ,
                f'overlap_mode "partial" and "containing" require non-empty anchor: {self}')

        assert self.fixed_placement is None or self.overlap_mode is None
        assert self.genotype and sum(self.genotype)
        for operation in self.operations:
            if operation.target_insertion_order is not None:
                operation.target_insertion_order = (self.sv_id,) + operation.target_insertion_order

    @property
    def breakends(self):
        return tuple(range(len(self.breakend_interval_lengths) + 1))

    @abstractmethod
    def to_vcf_records(self, config):
        raise NotImplementedError()

    def get_anchor_length(self):
        anchor = self.anchor
        if (self.overlap_mode == OverlapMode.EXACT) or (anchor is None):
            return None
        return sum(self.breakend_interval_lengths[anchor.start_breakend:anchor.end_breakend])

    def __str__(self):
        return self.config_descr

    def set_placement(self, placement, roi):
        placement = copy(placement)

        self.placement = placement
        self.roi = roi

        for operation in self.operations:
            operation.placement = placement

    def is_placed(self):
        return self.placement is not None

    def get_regions(self, placement=None):
        if placement is None:
            placement = self.placement

        regions = []
        for operation in self.operations:
            source_region = operation.get_source_region(placement)
            target_region = operation.get_target_region(placement)
            if source_region is not None:
                regions.append(source_region)
            if (target_region is not None) and (target_region != source_region):
                regions.append(target_region)
        regions_merged = []
        for region in sorted(set(regions)):
            if regions_merged and (region.chrom == regions_merged[-1].chrom) and (region.start == regions_merged[-1].end):
                regions_merged[-1] = regions_merged[-1].replace(end=region.end)
            else:
                regions_merged.append(region)

        return regions_merged

# end: class SV

class VariantType(Enum):
    INS = "INS"

    DEL = "DEL"
    INV = "INV"
    DUP = "DUP"
    mCNV = "mCNV"
    INV_DUP = "INV_DUP"
    DUP_INV = "DUP_INV"

    dDUP = "dDUP"
    INV_dDUP = "INV_dDUP"
    dDUP_INV = "dDUP_INV"
    INV_rTRA = "INV_rTRA"
    nrTRA = "nrTRA"
    rTRA = "rTRA"
    INV_nrTRA = "INV_nrTRA"

    delINV = "delINV"
    INVdel = "INVdel"
    dupINV = "dupINV"
    INVdup = "INVdup"

    INS_iDEL = "INS_iDEL"
    dDUP_iDEL = "dDUP_iDEL"

    dupINVdup = "dupINVdup"
    delINVdel = "delINVdel"
    delINVdup = "delINVdup"
    dupINVdel = "dupINVdel"

    SNP = "SNP"
    DIVERGENCE = "DIVERGENCE"

    Custom = "Custom"

    trEXP = "trEXP"
    trCON = "trCON"


class BaseSV(SV):

    @override
    def to_vcf_records(self, config):
        sv_type_str = self.info['SVTYPE']
        assert self.placement is not None

        sv_vcf_recs: list[dict] = []
        for op_idx, operation in enumerate(self.operations):

            dispersion_target = None

            if (not operation.is_in_place and
                    operation.source_breakend_region is not None):
                # Identify operations giving the target of a dispersion
                dispersion_target = operation.target_region
                assert dispersion_target.is_empty()

            op_type_str = operation.transform_type.value
            if (operation.transform_type == TransformType.IDENTITY) and operation.is_in_place:
                if (operation.transform.divergence_prob == 1) and (self.breakend_interval_lengths[0] == 1):
                    op_type_str = 'SNP'
                else:
                    op_type_str = 'DIVERGENCE'

            sv_info = dict(self.info)
            sv_info['IN_PLACE'] = "in_place" if operation.is_in_place else "paste"
            if operation.novel_insertion_seq is not None:
                op_chrom = operation.target_region.chrom
                op_start = operation.target_region.start
                op_end = op_start + 1
                svlen = len(operation.novel_insertion_seq)
                if config.get('output_vcf_ins_seq', True):
                    sv_info['INSSEQ'] = operation.novel_insertion_seq
            elif operation.source_region is not None:
                op_chrom = operation.source_region.chrom
                op_start = operation.source_region.start
                svlen = operation.source_region.end - operation.source_region.start
                op_end = operation.source_region.start + svlen
            else:
                op_chrom = operation.target_region.chrom
                op_start = operation.target_region.start
                op_end = op_start + 1
                svlen = 0

            if operation.transform.n_copies > 1:
                sv_info['NCOPIES'] = operation.transform.n_copies

            sv_info['SVTYPE'] = op_type_str
            if dispersion_target is not None:
                sv_info['TARGET_CHROM'] = dispersion_target.chrom
                sv_info['TARGET'] = dispersion_target.start + 1
            sv_info['SVLEN'] = svlen

            if operation.target_insertion_order is not None:
                sv_info['INSORD'] = operation.target_insertion_order[1]

            if self.anchor is not None:
                # If the operation breakend are within the anchor we give the overlap information.
                if (operation.source_breakend_region is not None and
                        self.anchor.start_breakend <= operation.source_breakend_region.start_breakend <= self.anchor.end_breakend):
                    sv_info['OVLP'] = self.roi.kind
                # Here the target is in the anchor.
                if (operation.target_region is not None and
                        self.anchor.start_breakend <= operation.target_region.start <= self.anchor.end_breakend):
                    sv_info['OVLP_TARGET'] = self.roi.kind
            sv_id = self.sv_id
            alleles = ['N', '<%s>' % op_type_str]
            if len(self.operations) > 1:
                sv_id += f'_{op_idx}'
                sv_info['PARENT_SVID'] = self.sv_id
                sv_info['PARENT_SVTYPE'] = sv_type_str
            else:
                sv_info['SVTYPE'] = sv_type_str
                alleles = ['N', '<%s>' % sv_type_str]
            for key, value in operation.op_info.items():
                sv_info[key] = value

            zygosity = tuple(map(int, self.genotype))
            vcf_rec = dict(contig=op_chrom, start=op_start, stop=op_end,
                           qual=100, filter='PASS',
                           alleles=alleles,
                           id=sv_id,
                           samples=[{'GT': zygosity}],
                           info=sv_info)
            sv_vcf_recs.append(vcf_rec)
        # end: for op_idx, operation in enumerate(self.operations)
        return sv_vcf_recs
# end: class BaseSV(SV)

class Syntax:
    DISPERSION = '_'

    DIVERGENCE = '*'
    MULTIPLE_COPIES = '+'

    ANCHOR_START = '('
    ANCHOR_END = ')'


SV_KEY = {
    VariantType.INS: ((), ("A",)),

    VariantType.DEL: (("A",), ()),
    VariantType.INV: (("A",), ("a",)),
    VariantType.DUP: (("A",), ("A", "A+")),
    VariantType.mCNV: (("A",), ("A+",)),
    VariantType.INV_DUP: (("A",), ("A", "a")),
    VariantType.DUP_INV: (("A",), ("a", "a")),

    VariantType.dDUP: (("A", "_"), ("A", "_", "A")),
    VariantType.INV_dDUP: (("A", "_"), ("A", "_", "a")),
    VariantType.dDUP_INV: (("A", "_"), ("a", "_", "a")),
    VariantType.INV_nrTRA: (("A", "_"), ("_", "a")),
    VariantType.nrTRA: (("A", "_"), ("_", "A")),
    VariantType.rTRA: (("A", "_", "B"), ("B", "_", "A")),
    VariantType.INV_rTRA: (("A", "_", "B"), ("b", "_", "a")),

    VariantType.delINV: (("A", "B"), ("b",)),
    VariantType.INVdel: (("A", "B"), ("a",)),
    VariantType.dupINV: (("A", "B"), ("A", "b", "a")),
    VariantType.INVdup: (("A", "B"), ("b", "a", "B")),

    VariantType.INS_iDEL: (("A", "_", "B"), ("_", "A")),
    VariantType.dDUP_iDEL: (("A", "_", "B"), ("A", "_", "A")),

    VariantType.dupINVdup: (("A", "B", "C"), ("A", "c", "b", "a", "C")),
    VariantType.delINVdel: (("A", "B", "C"), ("b",)),
    VariantType.delINVdup: (("A", "B", "C"), ("c", "b", "C")),
    VariantType.dupINVdel: (("A", "B", "C"), ("A", "b", "a")),

    VariantType.DIVERGENCE: (("A",), ("A*",)),
    VariantType.SNP: (("A",), ("A*",)),
}


@dataclass(order=True, frozen=True)
class Symbol:
    """A symbol denoting a region.

    Examples:

    A
    _1

    """

    name: str

    def __str__(self):
        return self.name


@dataclass(frozen=True)
class RHSItem:
    symbol: Symbol
    transform: Transform

    def __str__(self):
        rep = self.symbol.name
        if self.transform.transform_type == TransformType.INV:
            rep = rep.lower()
        if self.transform.divergence_prob > 0 or self.transform.replacement_seq:
            rep += Syntax.DIVERGENCE
        if self.transform.n_copies > 1:
            rep += Syntax.MULTIPLE_COPIES
        return rep

