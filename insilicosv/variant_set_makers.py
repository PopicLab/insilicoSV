from abc import ABC, abstractmethod
from collections import defaultdict
from contextlib import closing
import copy
from dataclasses import dataclass
import dataclasses
from enum import Enum
from itertools import chain, groupby
import math
import logging
import random
from typing_extensions import Optional, Any, ClassVar, Union, Type, Tuple, override

from pysam import FastaFile, VariantFile

from insilicosv import utils
from insilicosv.utils import (
    Region, RegionFilter, OverlapMode, Locus, error_context, chk,
    has_duplicates, pairwise)
from insilicosv.sv_defs import (Transform, TransformType,
                                Breakend, BreakendRegion, Operation, SV)

logger = logging.getLogger(__name__)

class VariantSetMaker(ABC):

    @classmethod
    @abstractmethod
    def can_make_from(cls, vset_config) -> bool:
        raise NotImplementedError()

    @abstractmethod
    def make_variant_set(self) -> list[SV]:
        raise NotImplementedError()

    @classmethod
    def get_vcf_header_infos(cls) -> list[dict]:
        return []

    next_sv_id: ClassVar[int] = 0

    @staticmethod
    def make_sv_id() -> str:
        sv_id = f'sv{VariantSetMaker.next_sv_id}'
        VariantSetMaker.next_sv_id += 1
        return sv_id

    def __init__(self, vset_config, sim_settings):
        self.vset_config = copy.deepcopy(vset_config)
        self.sim_settings = sim_settings

        self.preprocess_config()

    def preprocess_config(self) -> None:

        vset_cfg = self.vset_config
        vset_cfg['config_descr'] = str(vset_cfg)

        chk(isinstance(vset_cfg, dict), f'Each variant set must be a dict')

    def get_sampled_int_value(self, value, locals_dict: Optional[dict] = None) -> int:
        if isinstance(value, int):
            return value
        elif isinstance(value, list):
            chk(len(value) == 2, f'Expected [min, max] pair: {value}')
            chk(isinstance(value[0], int) and isinstance(value[1], int) and
                int(value[0]) <= int(value[1]), f'Invalid [min, max] pair: {value}')
            return random.randint(value[0], value[1])
        elif isinstance(value, str):
            try:
                eval_dict = dict(
                    random=random, math=math
                    )
                if locals_dict is not None:
                    eval_dict.update(locals_dict)
                return int(eval(value, eval_dict))
            except Exception as exc:
                chk(False, f'Error valuating {value}: {exc}')

class VariantType(Enum):
    INS = "INS"

    DEL = "DEL"
    INV = "INV"
    DUP = "DUP"
    INV_DUP = "INV_DUP"
    INV_DUP2 = "INV_DUP2"
    INV_DUP3 = "INV_DUP3"

    dDUP = "dDUP"
    INV_dDUP = "INV_dDUP"
    INV_dDUP2 = "INV_dDUP2"
    INV_dDUP3 = "INV_dDUP3"
    INV_TRA = "INV_TRA"
    TRA_NONRECIPROCAL = "TRA_NONRECIPROCAL"

    delINV = "delINV"
    INVdel = "INVdel"
    fldup_INV = "fldup_INV"
    INV_fldup = "INV_fldup"

    TRA_RECIPROCAL = "TRA_RECIPROCAL"
    INS_iDEL = "INS_iDEL"
    dDUP_iDEL = "dDUP_iDEL"

    dupINVdup = "dupINVdup"
    delINVdel = "delINVdel"
    delINVdup = "delINVdup"
    dupINVdel = "dupINVdel"

    SNP = "SNP"
    DIVERGENCE = "DIVERGENCE"

    Custom = "Custom"

    trINS = "trINS"
    trEXP = "trEXP"
    trCON = "trCON"

    TRA_CHROM_BALANCED = "TRA_CHROM_BALANCED"
    TRA_CHROM_UNBALANCED = "TRA_CHROM_UNBALANCED"

class BaseSV(SV):
    
    @override
    def to_vcf_records(self, sim_settings) -> list[dict]:
        sv_type_str = self.info['SVTYPE']

        assert self.placement is not None

        breakend_regions_inverted = [
            operation.source_breakend_region for operation in self.operations
            if (operation.is_in_place and
                operation.transform_type == TransformType.INV)]
        breakend_regions_deleted = [
            operation.source_breakend_region for operation in self.operations
            if operation.transform_type == TransformType.DEL]
        breakend_regions_changed = breakend_regions_inverted + breakend_regions_deleted
        breakend_regions_new_place = [
            operation.source_breakend_region for operation in self.operations
            if not operation.is_in_place]

        sv_vcf_recs: list[dict] = []
        for op_idx, operation in enumerate(self.operations):

            dispersion_target: Optional[Region] = None

            if (not operation.is_in_place and
                operation.source_breakend_region is not None and
                operation.target_insertion_breakend not in [
                    operation.source_breakend_region.start_breakend,
                    operation.source_breakend_region.end_breakend]):
                dispersion_target = operation.target_region
                assert dispersion_target.is_empty()

            op_type_str = ''

            if operation.is_in_place:
                if (operation.transform_type == TransformType.IDENTITY or
                    operation.source_breakend_region in breakend_regions_new_place):
                    continue
                if operation.transform_type == TransformType.DEL:
                    op_type_str = 'DEL'
                elif operation.transform_type == TransformType.INV:
                    op_type_str = 'INV'
            else:
                assert not operation.is_in_place
                if operation.novel_insertion_seq is not None:
                    op_type_str = 'INS'
                elif operation.chromosomal_translocation_source_breakend is not None:
                    if not operation.chromosomal_translocation_active:
                        continue
                    op_type_str = 'TRA_CHROMOSOMAL'
                else:
                    if dispersion_target is None:
                        if operation.source_breakend_region not in breakend_regions_changed:
                            if operation.transform_type == TransformType.IDENTITY:
                                op_type_str = 'DUP'
                            elif operation.transform_type == TransformType.INV:
                                if operation.target_insertion_breakend == operation.source_breakend_region.end_breakend:
                                    op_type_str = 'INV_DUP'  # A -> Aa'
                                elif operation.target_insertion_breakend == operation.source_breakend_region.start_breakend:
                                    op_type_str = 'INV_DUP2'  # A -> a'A
                                else:
                                    assert False
                        else:
                            assert operation.source_breakend_region in breakend_regions_changed
                            if operation.source_breakend_region in breakend_regions_inverted:
                                if operation.transform_type  == TransformType.INV:
                                    op_type_str = 'INV_DUP3'  # A -> aa'
                                else:
                                    chk(False, "A -> aA' should be written A -> a'A")
                            else:
                                assert operation.source_breakend_region in breakend_regions_deleted
                                chk(False, "A -> A' or A -> a' should not be needed")
                    else:
                        assert dispersion_target is not None
                        if operation.source_breakend_region not in breakend_regions_changed:
                            if operation.transform_type == TransformType.IDENTITY:
                                op_type_str = 'dDUP'
                            else:
                                op_type_str = 'INV_dDUP'
                        elif operation.source_breakend_region in breakend_regions_inverted:
                            if operation.transform_type == TransformType.IDENTITY:
                                op_type_str = 'INV_dDUP2'
                            else:
                                op_type_str = 'INV_dDUP3'
                        elif operation.source_breakend_region in breakend_regions_deleted:
                            if operation.transform_type == TransformType.IDENTITY:
                                op_type_str = 'TRA_NONRECIPROCAL'
                            else:
                                op_type_str = 'INV_TRA'

            assert op_type_str, (
                f'{operation=} {sv=} {dispersion_target=} '
                f'{breakend_regions_deleted=}')

            sv_info = dict(self.info)

            if operation.novel_insertion_seq is not None:
                op_chrom = operation.target_region.chrom
                op_start = operation.target_region.start
                op_end = op_start + 1
                svlen = len(operation.novel_insertion_seq)
                if sim_settings.get('output_vcf_ins_seq', True):
                    sv_info['INSSEQ'] = operation.novel_insertion_seq
            elif operation.chromosomal_translocation_source_breakend is not None:
                op_chrom = operation.target_region.chrom
                op_start = operation.target_region.start
                op_end = op_start + 1
                svlen = 0
            else:
                assert operation.source_region is not None
                op_chrom = operation.source_region.chrom
                op_start = operation.source_region.start
                op_end = operation.source_region.end
                svlen = op_end - op_start

            if operation.transform.n_copies > 1:
                sv_info['NCOPIES'] = operation.transform.n_copies

            sv_info['SVTYPE'] = op_type_str
            if dispersion_target is not None:
                sv_info['TARGET_CHROM'] = dispersion_target.chrom
                sv_info['TARGET'] = dispersion_target.start + 1
            sv_info['SVLEN'] = svlen

            if operation.target_insertion_order is not None:
                sv_info['INSORD'] = operation.target_insertion_order[1]

            if self.roi is not None and self.roi.kind and self.roi.kind != '_reference_':
                sv_info['OVLP'] = self.roi.kind

            zygosity = tuple(map(int, self.genotype))

            vcf_rec = dict(contig=op_chrom, start=op_start, stop=op_end,
                           qual=100, filter='PASS',
                           alleles=['N', '<%s>' % op_type_str],
                           id=self.sv_id,
                           samples=[{'GT': zygosity}],
                           info=sv_info)
            sv_vcf_recs.append(vcf_rec)
        # end: for op_idx, operation in enumerate(self.operations)

        if len(sv_vcf_recs) > 1:
            assert sv_type_str not in (
                'DEL', 'INV', 'INS', 'DUP', 'INV_DUP', 'INV_DUP2', 'INV_DUP3', 'dDUP', 'INV_dDUP',
                'INV_dDUP2', 'INV_dDUP3', 'TRA_NONRECIPROCAL', 'INV_TRA')
            for rec_num, sv_vcf_rec in enumerate(sv_vcf_recs):
                sv_vcf_rec['id'] += f'_{rec_num}'
                sv_vcf_rec['info']['PARENT_SVID'] = self.sv_id
                sv_vcf_rec['info']['PARENT_SVTYPE'] = sv_type_str

        return sv_vcf_recs
# end: class BaseSV(SV)

class SimulatedVariantSetMaker(VariantSetMaker):

    def __init__(self, vset_config, sim_settings):
        super().__init__(vset_config, sim_settings)

    @override
    def preprocess_config(self) -> None:
        super().preprocess_config()

        vset_cfg = self.vset_config
        chk(isinstance(vset_cfg.get('type'), str), f'Missing or bad variant type')
        chk(isinstance(vset_cfg.get('number'), int) and vset_cfg['number'] >= 0,
            f'Missing or bad "number" of variants to generate')

        try:
            vset_cfg['type'] = VariantType(vset_cfg['type'])
        except ValueError:
            chk(False, 'Invalid variant type')

        if 'overlap_mode' in vset_cfg:
            chk(isinstance(vset_cfg['overlap_mode'], str), 'overlap_mode not a string')
            try:
                vset_cfg['overlap_mode'] = OverlapMode(vset_cfg['overlap_mode'])
            except ValueError:
                chk(False, 'Invalid overlap_mode')

        if 'overlap_region_type' in vset_cfg:
            vset_cfg['overlap_region_type'] = utils.as_list(vset_cfg['overlap_region_type'])
        chk('overlap_region_type' not in vset_cfg or
            utils.is_list_of(str, vset_cfg['overlap_region_type']),
            'Invalid overlap_region_type')

        chk('divergence_prob' not in vset_cfg or
            (isinstance(vset_cfg['divergence_prob'], (int, float)) and
             0 <= vset_cfg['divergence_prob'] <= 1),
            'Invalid divergence_prob')

    # end: def preprocess_config(self)
        
    @property
    def sv_type(self) -> VariantType:
        return self.vset_config['type']

    def get_roi_filter(self) -> Optional[RegionFilter]:

        vset_config = self.vset_config
        if ('overlap_region_type' not in vset_config and
            'overlap_region_length_range' not in vset_config):
            return None

        return RegionFilter(
            region_kinds=(tuple(utils.as_list(vset_config['overlap_region_type']))
                          if 'overlap_region_type' in vset_config else None),
            region_length_range=
            tuple(vset_config.get('overlap_region_length_range', (None, None))))

    def get_blacklist_filter(self) -> Optional[RegionFilter]:

        vset_config = self.vset_config
        if 'blacklist_region_type' not in vset_config:
            return None

        return RegionFilter(
            region_kinds=tuple(utils.as_list(self.vset_config['blacklist_region_type'])))


    def pick_genotype(self) -> tuple[bool, bool]:
        if self.sim_settings.get('homozygous_only', False) or random.randint(0, 1):
            return (True, True)
        else:
            return random.choice([(True, False), (False, True)])

    @override
    def make_variant_set(self) -> list[SV]:
        return [self.simulate_sv() for _ in range(self.vset_config['number'])]

    @abstractmethod
    def simulate_sv(self) -> SV:
        raise NotImplementedError()

#############################################
# Constructing SVs from grammar definitions #
#############################################

class Syntax:
    DISPERSION = '_'

    DIVERGENCE = '*'
    MULTIPLE_COPIES = '+'
    NEW_PLACE = "^"

    ANCHOR_START = '('
    ANCHOR_END = ')'

    # for internal use only
    END_SYMBOL = '~'

SV_KEY: dict[VariantType, tuple[tuple[str, ...], tuple[str, ...]]] = {
    VariantType.INS: ((), ("A",)),

    VariantType.DEL: (("A",), ()),
    VariantType.INV: (("A",), ("a",)),
    VariantType.DUP: (("A",), ("A", "A+")),
    VariantType.INV_DUP: (("A",), ("A", "a^")),
    VariantType.INV_DUP2: (("A",), ("a^", "A")),
    VariantType.INV_DUP3: (("A",), ("a", "a^")),

    VariantType.dDUP: (("A", "_"), ("A", "_", "A^")),
    VariantType.INV_dDUP: (("A", "_"), ("A", "_", "a^")),
    VariantType.INV_dDUP2: (("A", "_"), ("a", "_", "A^")),
    VariantType.INV_dDUP3: (("A", "_"), ("a", "_", "a^")),
    VariantType.INV_TRA: (("A", "_"), ("_", "a^")),
    VariantType.TRA_NONRECIPROCAL: (("A", "_"), ("_", "A^")),

    VariantType.delINV: (("A", "B"), ("b",)),
    VariantType.INVdel: (("A", "B"), ("a",)),
    VariantType.fldup_INV: (("A", "B"), ("A", "b", "a^")),
    VariantType.INV_fldup: (("A", "B"), ("b^", "a", "B")),

    VariantType.TRA_RECIPROCAL: (("A", "_", "B"), ("B^", "_", "A^")),
    VariantType.INS_iDEL: (("A", "_", "B"), ("_", "A^")),
    VariantType.dDUP_iDEL: (("A", "_", "B"), ("A", "_", "A^")),

    VariantType.dupINVdup: (("A", "B", "C"), ("A", "c^", "b", "a^", "C")),
    VariantType.delINVdel: (("A", "B", "C"), ("b",)),
    VariantType.delINVdup: (("A", "B", "C"), ("c^", "b", "C")),
    VariantType.dupINVdel: (("A", "B", "C"), ("A", "b", "a^")),

    VariantType.DIVERGENCE: (("A",), ("A*",)),
    VariantType.SNP: (("A",), ("A*",)),
}

def n_dispersions(strs):
    return strs.count(Syntax.DISPERSION)

@dataclass(order=True, frozen=True)
class Symbol:

    """A symbol denoting a region.

    Examples: 

    A
    _1

    """
    
    name: str

    # END_SYMBOL: A special symbol representing empty sequence,
    # added to the end of LHS and RHS of each grammar rule to
    # avoid corner cases during processing.
    END_SYMBOL: ClassVar['Symbol'] = None  # type: ignore

    def is_regular(self) -> bool:
        return len(self.name) == 1 and self.name[0].isupper()

    def is_dispersion(self) -> bool:
        return (len(self.name) == 2 and self.name[0] == Syntax.DISPERSION
                and self.name[1].isdigit())

    @property
    def dispersion_number(self) -> int:
        assert self.is_dispersion()
        return int(self.name[1])

    def is_end_symbol(self) -> bool:
        return self.name == Syntax.END_SYMBOL

    def __post_init__(self) -> None:
        assert self.is_regular() or self.is_dispersion() or self.is_end_symbol()

    def __str__(self):
        return self.name

Symbol.END_SYMBOL = Symbol(Syntax.END_SYMBOL)

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
        if not self.transform.is_in_place:
            rep += Syntax.NEW_PLACE
        return rep

class FromGrammarVariantSetMaker(SimulatedVariantSetMaker):
    """Constructs an SV from a grammar definition"""

    @override
    @classmethod
    def can_make_from(cls, vset_config) -> bool:
        return vset_config.get("type") in [k.value for k in SV_KEY.keys()] + ['Custom']

    @override
    def preprocess_config(self) -> None:
        super().preprocess_config()

        for vset_config_key in self.vset_config:
            chk(vset_config_key in (
                'type', 'number',
                'source', 'target',
                'length_ranges',
                'overlap_region_type', 'overlap_region_length_range',
                'overlap_mode',
                'blacklist_region_type',
                'divergence_prob',
                'num_copies',
                'interchromosomal',
                'config_descr', 
            ), f'invalid SV config key {vset_config_key}')

        vset_cfg = self.vset_config

        if vset_cfg['type'] not in (VariantType.SNP,):
            chk('length_ranges' in vset_cfg, 'Please specify length ranges')
            chk(isinstance(vset_cfg['length_ranges'], list), 'length_ranges must be a list')
            for length_range in vset_cfg['length_ranges']:
                chk(isinstance(length_range, str) or 
                    (isinstance(length_range, list) and len(length_range) == 2 and
                     isinstance(length_range[0], (type(None), int, str)) and
                     isinstance(length_range[1], (type(None), int, str))),
                    'invalid length_ranges')
                
        if vset_cfg['type'] == VariantType.SNP:
            chk(vset_cfg.get('length_ranges') in (None, [[1, 1]]),
                'length range for SNP can only be [1, 1]')
            chk('divergence_prob' not in vset_cfg or vset_cfg['divergence_prob'] == 1,
                'divergence prob for SNP can only be 1')
            vset_cfg['length_ranges'] = [[1, 1]]
            vset_cfg['divergence_prob'] = 1.0

        chk('overlap_region_type' not in self.vset_config or
            'overlap_mode' in self.vset_config,
            'overlap_region_type is specified but not overlap_mode')

        if 'overlap_mode' in self.vset_config and 'overlap_region_type' not in self.vset_config:
            self.vset_config['overlap_region_type'] = ['all']

    @property
    def is_interchromosomal(self) -> bool:
        return self.vset_config.get('interchromosomal', False)

    def get_lhs_rhs_strs(self) -> tuple[tuple[str, ...], tuple[str, ...]]:
        for src_trg in ('source', 'target'):
            chk(isinstance(self.vset_config.get(src_trg), (type(None), list, str)),
                f'{src_trg} must be a string or list of strings')
            chk(not isinstance(self.vset_config.get(src_trg), list) or
                utils.is_list_of(str, self.vset_config.get(src_trg)),
                f'elements of {src_trg} must be strings')
                 
        if isinstance(self.vset_config.get('source'), str):
            self.vset_config['source'] = list(self.vset_config['source'])
            
        if isinstance(self.vset_config.get('target'), str):
            rhs_strs_list: list[str] = []
            for c in self.vset_config['target']:
                if c in (Syntax.DIVERGENCE, Syntax.MULTIPLE_COPIES, Syntax.NEW_PLACE):
                    chk(rhs_strs_list and rhs_strs_list[-1] and rhs_strs_list[-1][0].isalpha(),
                        f'{c} must modify a symbol')
                    rhs_strs_list[-1] += c
                else:
                    rhs_strs_list.append(c)

            self.vset_config['target'] = rhs_strs_list

        if self.sv_type == VariantType.Custom:
            chk('source' in self.vset_config and 'target' in self.vset_config,
                'Please provide "source" and "target" grammar for custom SV')
            lhs_strs, rhs_strs = (
                tuple(self.vset_config['source']), tuple(self.vset_config['target']))
        else:
            lhs_strs, rhs_strs = SV_KEY[VariantType(self.sv_type)]
            if 'source' in self.vset_config:
                lhs_strs = tuple(self.vset_config['source'])
                chk(tuple(s for s in lhs_strs if s not in (Syntax.ANCHOR_START, Syntax.ANCHOR_END)) ==
                    SV_KEY[VariantType(self.sv_type)][0],
                    f'Source spec does not match built-in spec')
            if 'target' in self.vset_config:
                rhs_strs = tuple(self.vset_config['target'])
                def cleanup_rhs_str(rhs_str):
                    return (rhs_str
                            .replace(Syntax.DIVERGENCE, ''))
                chk(tuple(cleanup_rhs_str(s)
                          for s in rhs_strs if s not in (Syntax.ANCHOR_START, Syntax.ANCHOR_END)) ==
                    tuple(cleanup_rhs_str(s)
                          for s in SV_KEY[VariantType(self.sv_type)][1]),
                    f'Target spec does not match built-in spec')

        all_strs = lhs_strs + rhs_strs
        chk(all(isinstance(s, str) for s in all_strs), 'Invalid LHS or RHS spec')
        chk((Syntax.ANCHOR_START not in all_strs and Syntax.ANCHOR_END not in all_strs) or
            (all_strs.count(Syntax.ANCHOR_START) == all_strs.count(Syntax.ANCHOR_END) == 1 and
             all_strs.index(Syntax.ANCHOR_START) <= all_strs.index(Syntax.ANCHOR_END) and
             (Syntax.ANCHOR_START in rhs_strs or Syntax.ANCHOR_END in lhs_strs)),
            'Invalid anchor spec')

        lhs_strs = self.equalize_lhs_rhs_dispersions(lhs_strs, rhs_strs)
        assert n_dispersions(lhs_strs) == n_dispersions(rhs_strs)

        if (self.sv_type != VariantType.Custom and
            rhs_strs.count(Syntax.DISPERSION) == 1 and
            not self.is_interchromosomal and
            random.randint(0, 1)):
            def flip_anchor(val: str) -> str:
                if val == Syntax.ANCHOR_START:
                    return Syntax.ANCHOR_END
                if val == Syntax.ANCHOR_END:
                    return Syntax.ANCHOR_START
                return val
            lhs_strs = tuple(map(flip_anchor, lhs_strs))[::-1]
            rhs_strs = tuple(map(flip_anchor, rhs_strs))[::-1]

        return lhs_strs, rhs_strs

    def equalize_lhs_rhs_dispersions(self, lhs_strs: tuple[str, ...],
                                     rhs_strs: tuple[str, ...]) -> tuple[str, ...]:
        # Example of equalizing dispersions:
        # ABC -> A_B'_C' becomes ABC__ -> A_B'_C'
        # ABC -> A'_B_C' becomes _ABC_ -> A'_B_C'
        if n_dispersions(lhs_strs) != n_dispersions(rhs_strs):
            chk(n_dispersions(lhs_strs) == 0,
                f'If LHS and RHS dispersion numbers differ then LHS must have zero dispersions')
            lhs_regular_symbols = [lhs_str for lhs_str in lhs_strs if lhs_str.isupper()]
            rhs_in_place_symbol_inds = [i for i, rhs_str in enumerate(rhs_strs)
                                        if (rhs_str[0].upper() in lhs_regular_symbols 
                                            and Syntax.NEW_PLACE in rhs_str)]
            chk(rhs_in_place_symbol_inds,
                f'RHS must have at least one in-place symbol to infer LHS dispersions')

            n_dispersions_prepend = n_dispersions(rhs_strs[:rhs_in_place_symbol_inds[0]])
            n_dispersions_append = n_dispersions(rhs_strs[rhs_in_place_symbol_inds[-1]+1:])

            chk(n_dispersions_prepend + n_dispersions_append == n_dispersions(rhs_strs),
                f'Cannot unambiguously infer how to balance LHS and RHS dispersions')
            lhs_strs = ((Syntax.DISPERSION,) * n_dispersions_prepend + lhs_strs +
                        (Syntax.DISPERSION,) * n_dispersions_append)
        return lhs_strs

    def initialize_grammar_rule(self) -> tuple[list[Symbol], list[RHSItem], Optional[BreakendRegion]]:

        lhs_strs, rhs_strs = self.get_lhs_rhs_strs()

        #
        # Parse the LHS strings into Symbols, and RHS strings into RHSItems.
        # Create unique Symbols for dispersions.
        #
        # If the overlap anchor was specified, determine the start and end
        # breakends of the anchor.
        #

        overlap_anchor_bounds: list[Optional[Breakend]] = [None, None]

        lhs: list[Symbol] = []
        
        n_dispersions_lhs = 0

        for lhs_str in lhs_strs + (Symbol.END_SYMBOL.name,):
            if lhs_str == Syntax.ANCHOR_START:
                overlap_anchor_bounds[0] = len(lhs)
                continue

            if lhs_str == Syntax.ANCHOR_END:
                overlap_anchor_bounds[1] = len(lhs)
                continue

            if lhs_str == Syntax.DISPERSION:
                n_dispersions_lhs += 1
                lhs_str += str(n_dispersions_lhs)
            else:
                chk((len(lhs_str) == 1 and lhs_str[0].isupper()) or
                    lhs_str == Symbol.END_SYMBOL.name,
                    f'Invalid LHS symbol {lhs_str}')

            lhs_item = Symbol(lhs_str)
            chk(lhs_item not in lhs, f'Duplicate LHS symbol {lhs_str}')
            lhs.append(lhs_item)
        
        rhs: list[RHSItem] = []
        n_dispersions_rhs = 0
        set_anchor_start = set_anchor_end = False
        last_in_place_rhs_item: Optional[RHSItem] = None
        for rhs_str in rhs_strs  + (Symbol.END_SYMBOL.name,):
            if rhs_str == Syntax.ANCHOR_START:
                set_anchor_start = True
                continue

            if rhs_str == Syntax.ANCHOR_END:
                set_anchor_end = True
                continue

            if rhs_str == Syntax.DISPERSION:
                n_dispersions_rhs += 1
                rhs_str += str(n_dispersions_rhs)
                symbol = Symbol(rhs_str)
            else:
                symbol = Symbol(rhs_str[0].upper())

            transform_type = TransformType.IDENTITY
            if rhs_str[0].islower():
                transform_type = TransformType.INV

            chk(Syntax.DIVERGENCE not in rhs_str or
                'divergence_prob' in self.vset_config,
                'Please specify divergence_prob')

            n_copies = 1
            if (Syntax.MULTIPLE_COPIES in rhs_str and
                'num_copies' in self.vset_config):
                n_copies = self.get_sampled_int_value(self.vset_config['num_copies'])

            is_in_place = (symbol in lhs and
                           Syntax.NEW_PLACE not in rhs_str and
                           Syntax.MULTIPLE_COPIES not in rhs_str and
                           sum(lhs_item.is_dispersion() for lhs_item in lhs[:lhs.index(symbol)]) ==
                           sum(rhs_item.symbol.is_dispersion() for rhs_item in rhs))

            transform = Transform(
                transform_type=transform_type,
                is_in_place=is_in_place,
                divergence_prob=(0 if Syntax.DIVERGENCE not in rhs_str else 
                                 self.vset_config.get('divergence_prob', 0)),
                n_copies=n_copies
            )
            rhs_item = RHSItem(symbol=symbol, transform=transform)
            rhs.append(rhs_item)

            if transform.is_in_place:
                chk(last_in_place_rhs_item is None or
                    lhs.index(symbol) > lhs.index(last_in_place_rhs_item.symbol),
                    f'RHS in-place item order must follow LHS order, but '
                    f'{rhs_item} comes after {last_in_place_rhs_item}')

                last_in_place_rhs_item = rhs_item

                if set_anchor_start:
                    overlap_anchor_bounds[0] = lhs.index(symbol)
                    set_anchor_start = False
                if set_anchor_end:
                    overlap_anchor_bounds[1] = lhs.index(symbol)
                    set_anchor_end = False

        overlap_anchor: Optional[BreakendRegion] = None
        if overlap_anchor_bounds[0] is not None and overlap_anchor_bounds[1] is not None:
            overlap_anchor = BreakendRegion(overlap_anchor_bounds[0], overlap_anchor_bounds[1])

        # If overlap_anchor is not specified, infer "full SV" as the anchor if possible
        if (overlap_anchor is None and
            'overlap_mode' in self.vset_config and not self.has_dispersions(lhs) and
            (self.vset_config['overlap_mode'] != OverlapMode.EXACT or
             (len(lhs) == 2 and lhs[0].is_regular())) and
            (len(lhs) > 1 or self.vset_config['overlap_mode'] == OverlapMode.CONTAINED)):
            overlap_anchor = BreakendRegion(0, len(lhs)-1)

        assert lhs and rhs and rhs[-1].transform.is_in_place
        assert not has_duplicates(lhs)
        rhs_in_place = [rhs_item.symbol for rhs_item in rhs if rhs_item.transform.is_in_place]
        assert not has_duplicates(rhs_in_place)
        assert all(lhs.index(s1) < lhs.index(s2)
                   for s1, s2 in utils.pairwise(rhs_in_place))

        return (lhs, rhs, overlap_anchor)

    # end: def initialize_grammar_rule(self)

    def novel_insertion_symbols(self, lhs: list[Symbol], rhs: list[RHSItem]) -> set[Symbol]:
        rhs_symbols = {rhs_item.symbol for rhs_item in rhs}
        return set(rhs_symbols) - set(lhs)

    def has_dispersions(self, lhs) -> bool:
        return any(symbol.is_dispersion() for symbol in lhs)

    def pick_symbol_lengths(
            self, lhs: list[Symbol], rhs: list[RHSItem]) -> Tuple[dict[Symbol, int], dict[Symbol, int]]:
        symbols_ordered = sorted(
            (symbol for symbol in chain(lhs, self.novel_insertion_symbols(lhs, rhs))
             if symbol.is_regular() or symbol.is_dispersion()),
            key=lambda symbol: (symbol.is_dispersion(),
                                symbol.name if not symbol.is_dispersion() else symbol.dispersion_number))

        chk('length_ranges' in self.vset_config, f'Need to specify length_ranges')
        length_ranges = self.vset_config['length_ranges']
        chk(isinstance(length_ranges, list) and len(length_ranges) == len(symbols_ordered),
            f'Expected {len(symbols_ordered)} length_ranges')
        chk(all(isinstance(length_range, str) or
                (isinstance(length_range, list) and len(length_range) == 2) for length_range in length_ranges),
            f'length_ranges must be a list of [min,max] pairs or expressions')

        symbol_lengths: dict[Symbol, int] = {Symbol.END_SYMBOL: 0}
        symbol_min_lengths: dict[Symbol, int] = {}
        for compute_dependent_lengths in (False, True):
            for symbol, length_range in zip(symbols_ordered, length_ranges):
                if symbol in symbol_lengths:
                    assert compute_dependent_lengths
                    continue

                if isinstance(length_range, (int, str)):
                    length_val = self.get_sampled_int_value(length_range)
                    length_range = (length_val, length_val)

                min_length, max_length = length_range

                if isinstance(min_length, str) or isinstance(max_length, str):
                    if not compute_dependent_lengths:
                        continue
                    expr_namespace: dict[str, int] = {
                        symbol.name: symbol_length for symbol, symbol_length in symbol_lengths.items()
                        if symbol.is_regular()}
                    try:
                        min_length = int(eval(str(min_length), expr_namespace))
                        max_length = int(eval(str(max_length), expr_namespace))
                    except Exception as exc:
                        chk(False, f'Error evaluating dependent length expr: {exc}')
                chk(isinstance(min_length, (int, type(None))) and
                    isinstance(max_length, (int, type(None))),
                    f'min, max must be int or null')

                if max_length is None:
                    if min_length is not None:
                        chk(symbol.is_dispersion(), f'Only dispersions can have min but not max length')
                        symbol_min_lengths[symbol] = min_length
                else:
                    chk(min_length is not None, f'max_length given but not min_length')
                    chk(min_length <= max_length,
                        f'max length less than min length')
                    chk(min_length >= 0, 'min length cannot be negative')
                    symbol_lengths[symbol] = random.randint(min_length, max_length)

        return symbol_lengths, symbol_min_lengths

    @override
    def simulate_sv(self) -> SV:

        lhs: list[Symbol]
        rhs: list[RHSItem]
        anchor: Optional[BreakendRegion]

        lhs, rhs, anchor = self.initialize_grammar_rule()

        symbol_lengths: dict[Symbol, int]
        symbol_min_lengths: dict[Symbol, int]
        symbol_lengths, symbol_min_lengths = self.pick_symbol_lengths(lhs, rhs)

        chk(not self.is_interchromosomal or
            (self.has_dispersions(lhs) and
             all(symbol not in symbol_lengths and symbol not in symbol_min_lengths
                 for symbol in lhs if symbol.is_dispersion())),
            f'Interchromosomal dispersions must be unbounded')

        novel_insertions: dict[Symbol, str] = {  # type: ignore
            symbol: utils.generate_seq(length=symbol_lengths[symbol])
            for symbol in self.novel_insertion_symbols(lhs, rhs)}

        breakend_interval_lengths = [symbol_lengths.get(lhs_item) for lhs_item in lhs]
        breakend_interval_min_lengths = [symbol_min_lengths.get(lhs_item) for lhs_item in lhs]
        symbol_start_breakends: dict[Symbol, Breakend] = {
            lhs_item: lhs_item_index for lhs_item_index, lhs_item in enumerate(lhs)}
        operations: list[Operation] = []

        assert lhs and rhs and rhs[-1].transform.is_in_place

        for rhs_item_idx in range(len(rhs)-1, -1, -1):
            rhs_item = rhs[rhs_item_idx]
            if rhs_item.transform.is_in_place:
                current_breakend = symbol_start_breakends[rhs_item.symbol]

            if rhs_item.symbol.is_regular():
                operation = Operation(transform=rhs_item.transform)
                if rhs_item.symbol in self.novel_insertion_symbols(lhs, rhs):
                    operation.novel_insertion_seq = novel_insertions[rhs_item.symbol]
                else:
                    operation.source_breakend_region = BreakendRegion(
                        symbol_start_breakends[rhs_item.symbol],
                        symbol_start_breakends[rhs_item.symbol]+1)

                if not rhs_item.transform.is_in_place:
                    operation.target_insertion_breakend = current_breakend
                    operation.target_insertion_order = (rhs_item_idx,)

                operation.op_info = dict(sv_type=self.sv_type,
                                         grammar_item=rhs_item,
                                         symbol=rhs_item.symbol)

                operations.append(operation)
        # end: for rhs_item_idx in range(len(rhs)-1, -1, -1)

        # add deletion operations for LHS symbols that have no in-place transform
        rhs_in_place_symbols = {
            rhs_item.symbol for rhs_item in rhs if rhs_item.transform.is_in_place}
        for lhs_item in lhs:
            if lhs_item not in rhs_in_place_symbols:
                operations.append(Operation(
                    source_breakend_region=BreakendRegion(symbol_start_breakends[lhs_item],
                                                          symbol_start_breakends[lhs_item]+1),
                    transform=Transform(transform_type=TransformType.DEL,
                                        is_in_place=True),
                    op_info=dict(sv_type=self.sv_type,
                                 grammar_item=lhs_item,
                                 symbol=lhs_item)))

        #
        # construct the SV object
        #

        info = self.construct_info()

        overlap_mode = self.vset_config.get('overlap_mode')
        roi_filter = self.get_roi_filter()

        if overlap_mode is None:
            assert roi_filter is None and anchor is None

            # unconstrained placement
            overlap_mode = OverlapMode.CONTAINED
            roi_filter = RegionFilter(region_kinds=('_reference_',))

            # find longest possible dispersion-free anchor
            anchors = [BreakendRegion(len(breakend_interval_lengths), len(breakend_interval_lengths))]
            for k, g in groupby(range(len(breakend_interval_lengths)),
                                key=lambda breakend: breakend_interval_lengths[breakend] is not None):
                if k:
                    breakends = list(g)
                    anchors.append(BreakendRegion(breakends[0], breakends[-1]+1))
            anchor = sorted(anchors, key=lambda breakend_region: sum(
                breakend_interval_lengths[
                    breakend_region.start_breakend:breakend_region.end_breakend]))[-1]  # type: ignore


        chk(overlap_mode is None or anchor is not None,
            f'overlap_mode specified but anchor neither specified nor inferrable')


        return BaseSV(sv_id=self.make_sv_id(),
                      breakend_interval_lengths=breakend_interval_lengths,
                      breakend_interval_min_lengths=breakend_interval_min_lengths,
                      is_interchromosomal=self.is_interchromosomal,
                      operations=operations,
                      anchor=anchor,
                      overlap_mode=overlap_mode,
                      roi_filter=roi_filter,
                      blacklist_filter=self.get_blacklist_filter(),
                      fixed_placement=None,
                      info=info,
                      genotype=self.pick_genotype(),
                      config_descr=self.vset_config['config_descr'])

    def construct_info(self) -> dict[str, Any]:
        sv_type_str = self.sv_type.value
        if self.sv_type == VariantType.Custom:
            source_str = ''.join(self.vset_config['source'])
            target_str = ''.join(self.vset_config['target'])
            sv_type_str += f'_{source_str}>{target_str}'
        return {'SVTYPE': sv_type_str}

# end: class FromGrammarVariantSetMaker

#############################################
# Constructing SVs involving tandem repeats #
#############################################

@dataclass(frozen=True)
class TandemRepeatRegionFilter(RegionFilter):
    min_num_repeats: int = 0

    @override
    def satisfied_for(self, region: Region) -> bool:
        if not super().satisfied_for(region):
            return False
        try:
            repeat_unit_length: int = int(region.data)
        except ValueError:
            return False

        return repeat_unit_length * self.min_num_repeats <= region.length()

@dataclass
class TandemRepeatExpansionContractionSV(BaseSV):

    num_repeats_in_placement: int = 0

    @override
    def set_placement(self, placement: list[Locus], roi: Optional[Region]) -> None:
        assert len(placement) == 2
        assert roi is not None
        repeat_unit_length: int = int(roi.data)
        assert 0 < repeat_unit_length <= roi.length()
        shift_length: int = repeat_unit_length * self.num_repeats_in_placement
        assert shift_length <= roi.length()
        super().set_placement(placement=[placement[0], placement[0].shifted(shift_length)],
                              roi=roi)

    @override
    def to_vcf_records(self, sim_settings) -> list[dict]:
        sv_type_str = self.info['SVTYPE']
        vcf_rec = super().to_vcf_records(sim_settings)[0]
        vcf_rec['info']['SVTYPE'] = sv_type_str
        vcf_rec['alleles'] = ['N', '<%s>' % sv_type_str]
        assert self.roi is not None
        vcf_rec['info']['SVLEN'] = self.roi.length()
        vcf_rec['stop'] = self.roi.end
        return [vcf_rec]
                
class TandemRepeatVariantSetMaker(SimulatedVariantSetMaker):

    @override
    @classmethod
    def can_make_from(cls, vset_config) -> bool:
        return vset_config.get("type") in ('trINS', 'trEXP', 'trCON')

    def __init__(self, vset_config, sim_settings):
        super().__init__(vset_config, sim_settings)

        chk('repeat_count_change_range' in self.vset_config,
            'Please specify repeat_count_change_range')
        chk(utils.is_valid_closed_int_range(self.vset_config['repeat_count_change_range']),
            'Invalid repeat_count_change_range')

        if self.vset_config["type"] == VariantType.trINS:
            chk('repeat_sequence_choices' in self.vset_config,
                'Please specify repeat_sequence_choices')
            chk(self.vset_config['repeat_sequence_choices'] and
                isinstance(self.vset_config['repeat_sequence_choices'], list) and
                all(isinstance(s, str) and len(s) >= 1 and all(c in 'TCGAtcga' for c in s)
                    for s in self.vset_config['repeat_sequence_choices']),
                'Invalid repeat_sequence_choices')

    @override
    def preprocess_config(self) -> None:
        super().preprocess_config()

        for vset_config_key in self.vset_config:
            chk(vset_config_key in (
                'type', 'number',
                'repeat_count_change_range',
                'repeat_sequence_choices',
                'overlap_region_type', 'overlap_region_length_range',
                'overlap_mode',
                'blacklist_region_type',
                'config_descr', 
            ), f'invalid SV config key {vset_config_key}')

        if self.vset_config['type'] in (VariantType.trEXP, VariantType.trCON):
            chk('overlap_region_type' in self.vset_config,
                'Please specify overlap_region_type corresponding to existing tandem repeats')

    @override
    def simulate_sv(self) -> SV:

        repeat_count_change = random.randint(*self.vset_config['repeat_count_change_range'])

        overlap_mode: Optional[OverlapMode] = self.vset_config.get('overlap_mode')

        anchor: Optional[BreakendRegion] = None
        roi_filter = self.get_roi_filter()
        breakend_interval_lengths: list[Optional[int]]

        if self.sv_type == VariantType.trEXP:
            chk(overlap_mode in [None, OverlapMode.EXACT],
                'overlap_mode must be "exact" for trEXP')
            breakend_interval_lengths = [None]
            operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY,
                                                        is_in_place=False,
                                                        n_copies=repeat_count_change),
                                    source_breakend_region=BreakendRegion(start_breakend=0, end_breakend=1),
                                    target_insertion_breakend=0,
                                    target_insertion_order=(0,))]
            anchor = BreakendRegion(0, 1)
            overlap_mode = OverlapMode.EXACT
            info = dict(SVTYPE=self.sv_type.value, TR_CHANGE=repeat_count_change)
            return TandemRepeatExpansionContractionSV(
                sv_id=self.make_sv_id(),
                breakend_interval_lengths=breakend_interval_lengths,
                breakend_interval_min_lengths=[None] * len(breakend_interval_lengths),
                is_interchromosomal=False,
                operations=operations,
                anchor=anchor,
                overlap_mode=overlap_mode,
                roi_filter=roi_filter,
                blacklist_filter=self.get_blacklist_filter(),
                fixed_placement=None,
                info=info,
                genotype=self.pick_genotype(),
                config_descr=self.vset_config['config_descr'],
                num_repeats_in_placement=1)
        elif self.sv_type == VariantType.trCON:
            breakend_interval_lengths = [None]
            operations = [Operation(transform=Transform(transform_type=TransformType.DEL,
                                                        is_in_place=True),
                                    source_breakend_region=BreakendRegion(0, 1))]
            anchor = BreakendRegion(0, 1)
            overlap_mode = OverlapMode.EXACT

            # ensure there are enough existing repeats to delete.
            overlap_region_type = (tuple(utils.as_list(self.vset_config['overlap_region_type']))
                                   if 'overlap_region_type' in self.vset_config else None)
            roi_filter = TandemRepeatRegionFilter(min_num_repeats=repeat_count_change+1,
                                                  region_kinds=overlap_region_type)
            info = dict(SVTYPE=self.sv_type.value, TR_CHANGE=repeat_count_change)
            return TandemRepeatExpansionContractionSV(
                sv_id=self.make_sv_id(),
                breakend_interval_lengths=breakend_interval_lengths,
                breakend_interval_min_lengths=[None] * len(breakend_interval_lengths),
                is_interchromosomal=False,
                operations=operations,
                anchor=anchor,
                overlap_mode=overlap_mode,
                roi_filter=roi_filter,
                blacklist_filter=self.get_blacklist_filter(),
                fixed_placement=None,
                info=info,
                genotype=self.pick_genotype(),
                config_descr=self.vset_config['config_descr'],
                num_repeats_in_placement=repeat_count_change)

        elif self.sv_type == VariantType.trINS:
            breakend_interval_lengths = []

            repeat_sequence = random.choice(self.vset_config['repeat_sequence_choices'])

            operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY,
                                                        is_in_place=False),
                                    novel_insertion_seq=repeat_sequence * repeat_count_change,
                                    target_insertion_breakend=0,
                                    target_insertion_order=(0,))]

            anchor = BreakendRegion(0, 0)
            if overlap_mode is None:
                assert roi_filter is None
                overlap_mode = OverlapMode.CONTAINED
                roi_filter = RegionFilter(region_kinds=('_reference_',))

            info = dict(SVTYPE=self.sv_type.value, TR_CHANGE=repeat_count_change,
                        TR_INS_REP_SEQ=repeat_sequence)

            return BaseSV(
                sv_id=self.make_sv_id(),
                breakend_interval_lengths=breakend_interval_lengths,
                breakend_interval_min_lengths=[None] * len(breakend_interval_lengths),
                is_interchromosomal=False,
                operations=operations,
                anchor=anchor,
                overlap_mode=overlap_mode,
                roi_filter=roi_filter,
                blacklist_filter=self.get_blacklist_filter(),
                fixed_placement=None,
                info=info,
                genotype=self.pick_genotype(),
                config_descr=self.vset_config['config_descr'])
        else:
            assert False

    @override
    @classmethod
    def get_vcf_header_infos(cls) -> list[dict]:
        return [dict(id='TR_INS_REP_SEQ', number=1, type='String',
                            description="for tandem repeat insertions, the sequence of a single repeat"),
                dict(id='TR_CHANGE', number=1, type='Integer',
                     description="for tandem repeat insertions/expansions/deletions, number of "
                     "repeats inserted or removed")]


##############################
# Chromosomal translocations #
##############################

class ChromosomalTranslocationVariantSetMaker(SimulatedVariantSetMaker):
    
    @override
    @classmethod
    def can_make_from(cls, vset_config) -> bool:
        return vset_config.get("type") in ('TRA_CHROM_BALANCED', 'TRA_CHROM_UNBALANCED')

    @override
    def preprocess_config(self) -> None:
        super().preprocess_config()

        for vset_config_key in self.vset_config:
            chk(vset_config_key in (
                'type', 'number',
                'overlap_region_type', 'overlap_region_length_range',
                'overlap_mode',
                'blacklist_region_type',
                'config_descr', 
            ), f'invalid SV config key {vset_config_key}')

    @override
    def simulate_sv(self) -> SV:
        overlap_mode: Optional[OverlapMode] = self.vset_config.get('overlap_mode')
        anchor = BreakendRegion(0, 0)
        roi_filter = self.get_roi_filter()

        if overlap_mode is None:
            assert roi_filter is None
            overlap_mode = OverlapMode.CONTAINED
            roi_filter = RegionFilter(region_kinds=('_reference_',))


        operations: list[Operation] = [Operation(transform=Transform(transform_type=TransformType.IDENTITY,
                                                                     is_in_place=False),
                                                 chromosomal_translocation_source_breakend=1,
                                                 chromosomal_translocation_active=True,
                                                 target_insertion_breakend=0,
                                                 target_insertion_order=(0,)),
                                       Operation(transform=Transform(transform_type=TransformType.IDENTITY,
                                                                     is_in_place=False),
                                                 chromosomal_translocation_source_breakend=0,
                                                 chromosomal_translocation_active=self.vset_config['type'] == VariantType.TRA_CHROM_BALANCED,
                                                 target_insertion_breakend=1,
                                                 target_insertion_order=(0,))]

        info = dict(SVTYPE=self.sv_type.value)

        return BaseSV(sv_id=self.make_sv_id(),
                      breakend_interval_lengths=[None],
                      breakend_interval_min_lengths=[None],
                      is_interchromosomal=True,
                      operations=operations,
                      anchor=anchor,
                      overlap_mode=overlap_mode,
                      roi_filter=roi_filter,
                      blacklist_filter=self.get_blacklist_filter(),
                      fixed_placement=None,
                      info=info,
                      genotype=self.pick_genotype(),
                      config_descr=self.vset_config['config_descr'])

###########################
# Importing SVs from VCFs #
###########################

class ImportedVariantSetMaker(VariantSetMaker):

    @override
    @classmethod
    def can_make_from(cls, vset_config) -> bool:
        return 'import' in vset_config

    def __init__(self, vset_config, sim_settings):
        super().__init__(vset_config, sim_settings)

        chk(utils.is_readable_file(vset_config['import']),
            'vcf must name a readable file')
        chk(set(vset_config.keys()) <= {'import'}, 'invalid config key')

        with FastaFile(sim_settings['reference']) as reference:
            self.chrom_lengths = { chrom: chrom_length
                                   for chrom, chrom_length in zip(reference.references,
                                                                  reference.lengths) }

    @override
    def make_variant_set(self) -> list[SV]:
        svs = []
        with closing(VariantFile(self.vset_config['import'])) as vcf:
            for vcf_rec in vcf.fetch():
                with error_context(vcf_rec):
                    sv = self.import_sv_from_vcf_rec(vcf_rec, self.sim_settings)
                    if sv is not None:
                        svs.append(sv)
        return self.restore_parent_svs(svs)

    def import_sv_from_vcf_rec(self, vcf_rec, sim_settings) -> Optional[SV]:

        chk(vcf_rec.chrom in self.chrom_lengths, f'unknown contig {vcf_rec.chrom}')
        chk(len(vcf_rec.alleles) == 2, 'Can only import bi-allelic variants from VCF')
        chk(len(vcf_rec.samples) == 1, 'Can only import VCFs with one sample')

        sample = vcf_rec.samples[0]
        genotype: tuple[bool, bool]
        if 'AF' in vcf_rec.info:
            alt_allele_freq = vcf_rec.info['AF'][0]
            chk(isinstance(alt_allele_freq, float) and 0 <= alt_allele_freq <= 1,
                f'invalid AF value {alt_allele_freq}')
            chk(not sample.phased, 'cannot use allele frequency-based import with phased sample')
            genotype = (random.random() < alt_allele_freq, random.random() < alt_allele_freq)
            if genotype == (False, False):
                return None
        elif sample['GT'] == (1, 1):
            genotype = (True, True)
        else:
            if not sample.phased:
                genotype = random.choice([(True, False), (False, True)])
            else:
                chk(sample['GT'] in ((1, 0), (0, 1)), 'invalid genotype')
                genotype = (bool(sample['GT'][0]), bool(sample['GT'][1]))

        vcf_info = dict(vcf_rec.info)
        if set(''.join(vcf_rec.alleles).upper()) <= set('TCGA'):
            if len(vcf_rec.alleles[0]) == 1 and len(vcf_rec.alleles[1]) == 1:
                # SNP
                vcf_info['SVTYPE'] = 'SNP'

        chk('SVTYPE' in vcf_info, 'Currently need SVTYPE to import from vcf')
        sv_type_str = vcf_info['SVTYPE']
        can_import_types = (VariantType.INS, VariantType.DEL, VariantType.INV, VariantType.DUP, VariantType.INV_DUP3,
                            VariantType.INV_DUP, VariantType.INV_DUP2, VariantType.dDUP, VariantType.INV_dDUP,
                            VariantType.INV_dDUP2, VariantType.INV_dDUP3, VariantType.TRA_NONRECIPROCAL,
                            VariantType.INV_TRA, VariantType.SNP)
        chk(sv_type_str in {variant_type.name for variant_type in can_import_types},
            f'Currently can only import from vcf these SV types: {can_import_types}')

        sv_type = VariantType(sv_type_str)

        sv_id = self.make_sv_id()

        sv_len = abs(vcf_info.get('SVLEN', vcf_rec.stop - vcf_rec.start))
        sv_start = Locus(chrom=vcf_rec.chrom, pos=vcf_rec.start)
        sv_end = sv_start.shifted(sv_len)

        target: Optional[Locus] = None
        is_interchromosomal = False
        if Syntax.DISPERSION in SV_KEY[sv_type][1]:
            chk(isinstance(vcf_info.get('TARGET'), int), 'need TARGET to import this vcf record')
            target = Locus(chrom=vcf_info.get('TARGET_CHROM', vcf_rec.chrom), pos=vcf_info['TARGET']-1)
            if target.chrom != vcf_rec.chrom:
                is_interchromosomal = True
            else:
                chk(target.pos >= sv_end.pos or target.pos <= sv_start.pos,
                    "invalid dispersion target")

        insord = vcf_info.get('INSORD', 0)

        breakend_interval_lengths: list[Optional[Breakend]]
        fixed_placement: list[Locus]

        if sv_type == VariantType.INS:
            if 'INSSEQ' in vcf_info:
                chk(isinstance(vcf_info['INSSEQ'], str), 'invalid INSSEQ type in vcf')
                novel_insertion_seq = vcf_info['INSSEQ']
            else:
                novel_insertion_seq = utils.generate_seq(length=sv_len)
            fixed_placement = [sv_start]
            operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY, is_in_place=False),
                                    novel_insertion_seq=novel_insertion_seq,
                                    target_insertion_breakend=0, target_insertion_order=(insord,))]
        elif sv_type == VariantType.DEL:
            fixed_placement = [sv_start, sv_end]
            operations = [Operation(transform=Transform(transform_type=TransformType.DEL, is_in_place=True),
                                    source_breakend_region=BreakendRegion(0, 1)),]
        elif sv_type == VariantType.INV:
            fixed_placement = [sv_start, sv_end]
            operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=True),
                                    source_breakend_region=BreakendRegion(0, 1)),]
        elif sv_type == VariantType.DUP:
            fixed_placement = [sv_start, sv_end]
            operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY, is_in_place=False,
                                                        n_copies=vcf_info.get('NCOPIES', 1)),
                                    source_breakend_region=BreakendRegion(0, 1),
                                    target_insertion_breakend=1, target_insertion_order=(insord,))]
        elif sv_type == VariantType.INV_DUP3:
            # A -> aa'
            fixed_placement = [sv_start, sv_end]
            operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=True),
                                    source_breakend_region=BreakendRegion(0, 1)),
                          Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                    source_breakend_region=BreakendRegion(0, 1),
                                    target_insertion_breakend=1, target_insertion_order=(insord,))]
        elif sv_type == VariantType.INV_DUP:
            # A -> Aa'
            fixed_placement = [sv_start, sv_end]
            operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                    source_breakend_region=BreakendRegion(0, 1),
                                    target_insertion_breakend=1, target_insertion_order=(insord,))]
        elif sv_type == VariantType.INV_DUP2:
            # A -> a'A
            fixed_placement = [sv_start, sv_end]
            operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                    source_breakend_region=BreakendRegion(0, 1),
                                    target_insertion_breakend=0, target_insertion_order=(insord,))]
        elif sv_type == VariantType.dDUP:
            # A_ -> A_A'
            assert target is not None
            if target.chrom != vcf_rec.chrom or target.pos >= sv_end.pos:
                fixed_placement = [sv_start, sv_end, target]
                operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY, is_in_place=False),
                                        source_breakend_region=BreakendRegion(0, 1),
                                        target_insertion_breakend=2, target_insertion_order=(insord,))]
            else:
                # _A -> A'_A
                fixed_placement = [target, sv_start, sv_end]
                operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY, is_in_place=False),
                                        source_breakend_region=BreakendRegion(1, 2),
                                        target_insertion_breakend=0, target_insertion_order=(insord,))]

        elif sv_type == VariantType.INV_dDUP:
            # A_ -> A_a'
            assert target is not None
            if target.chrom != vcf_rec.chrom or target.pos >= sv_end.pos:
                fixed_placement = [sv_start, sv_end, target]
                operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                        source_breakend_region=BreakendRegion(0, 1),
                                        target_insertion_breakend=2, target_insertion_order=(insord,))]
            else:
                # _A -> a'_A
                fixed_placement = [target, sv_start, sv_end]
                operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                        source_breakend_region=BreakendRegion(1, 2),
                                        target_insertion_breakend=0, target_insertion_order=(insord,))]

        elif sv_type == VariantType.INV_dDUP2:
            # A_ -> a_A'
            assert target is not None
            if target.chrom != vcf_rec.chrom or target.pos >= sv_end.pos:
                fixed_placement = [sv_start, sv_end, target]
                operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=True),
                                        source_breakend_region=BreakendRegion(0, 1)),
                              Operation(transform=Transform(transform_type=TransformType.IDENTITY, is_in_place=False),
                                        source_breakend_region=BreakendRegion(0, 1),
                                        target_insertion_breakend=2, target_insertion_order=(insord,))]
            else:
                # _A -> A'_a
                fixed_placement = [target, sv_start, sv_end]
                operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=True),
                                        source_breakend_region=BreakendRegion(1, 2)),
                              Operation(transform=Transform(transform_type=TransformType.IDENTITY, is_in_place=False),
                                        source_breakend_region=BreakendRegion(1, 2),
                                        target_insertion_breakend=0, target_insertion_order=(insord,))]

        elif sv_type == VariantType.INV_dDUP3:
            # A_ -> a_a'
            assert target is not None
            if target.chrom != vcf_rec.chrom or target.pos >= sv_end.pos:
                fixed_placement = [sv_start, sv_end, target]
                operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=True),
                                        source_breakend_region=BreakendRegion(0, 1)),
                              Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                        source_breakend_region=BreakendRegion(0, 1),
                                        target_insertion_breakend=2, target_insertion_order=(insord,))]
            else:
                # _A -> a'_a
                fixed_placement = [target, sv_start, sv_end]
                operations = [Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=True),
                                        source_breakend_region=BreakendRegion(1, 2)),
                              Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                        source_breakend_region=BreakendRegion(1, 2),
                                        target_insertion_breakend=0, target_insertion_order=(insord,))]
        elif sv_type == VariantType.TRA_NONRECIPROCAL:
            # A_ -> _A'
            assert target is not None
            if target.chrom != vcf_rec.chrom or target.pos >= sv_end.pos:
                fixed_placement = [sv_start, sv_end, target]
                operations = [Operation(transform=Transform(transform_type=TransformType.DEL, is_in_place=True),
                                        source_breakend_region=BreakendRegion(0, 1)),
                              Operation(transform=Transform(transform_type=TransformType.IDENTITY, is_in_place=False),
                                        source_breakend_region=BreakendRegion(0, 1),
                                        target_insertion_breakend=2, target_insertion_order=(insord,))]
            else:
                # _A -> A'_
                fixed_placement = [target, sv_start, sv_end]
                operations = [Operation(transform=Transform(transform_type=TransformType.DEL, is_in_place=True),
                                        source_breakend_region=BreakendRegion(1, 2)),
                              Operation(transform=Transform(transform_type=TransformType.IDENTITY, is_in_place=False),
                                        source_breakend_region=BreakendRegion(1, 2),
                                        target_insertion_breakend=0, target_insertion_order=(insord,))]
        elif sv_type == VariantType.INV_TRA:
            # A_ -> _a'
            assert target is not None
            if target.chrom != vcf_rec.chrom or target.pos >= sv_end.pos:
                fixed_placement = [sv_start, sv_end, target]
                operations = [Operation(transform=Transform(transform_type=TransformType.DEL, is_in_place=True),
                                        source_breakend_region=BreakendRegion(0, 1)),
                              Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                        source_breakend_region=BreakendRegion(0, 1),
                                        target_insertion_breakend=2, target_insertion_order=(insord,))]
            else:
                # _A -> a'_
                fixed_placement = [target, sv_start, sv_end]
                operations = [Operation(transform=Transform(transform_type=TransformType.DEL, is_in_place=True),
                                        source_breakend_region=BreakendRegion(1, 2)),
                              Operation(transform=Transform(transform_type=TransformType.INV, is_in_place=False),
                                        source_breakend_region=BreakendRegion(1, 2),
                                        target_insertion_breakend=0, target_insertion_order=(insord,))]
        elif sv_type == VariantType.SNP:
            # A -> A*
            breakend_interval_lengths = [1]
            chk(len(vcf_rec.alleles[1]) == 1, 'invalid alt allele length for a SNP')
            operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY,
                                                        is_in_place=True,
                                                        replacement_seq=vcf_rec.alleles[1]),
                                    source_breakend_region=BreakendRegion(0, 1))]
            fixed_placement = [sv_start, sv_end]
        else:
            chk(False, f'cannot import from VCF variants of type {sv_type}')

        valid_info_keys = {vcf_header_info['id'] for vcf_header_info in get_vcf_header_infos()} - {'END', 'TARGET', 'TARGET_CHROM', 'NCOPIES'}
        info = {key: val for key, val in vcf_info.items() if key in valid_info_keys}

        breakend_interval_lengths = [
            locus2.pos - locus1.pos if locus1.chrom == locus2.chrom else None
            for locus1, locus2 in pairwise(fixed_placement)]

        return BaseSV(sv_id=sv_id,
                      breakend_interval_lengths=breakend_interval_lengths,
                      breakend_interval_min_lengths=[None] * len(breakend_interval_lengths),
                      is_interchromosomal=False,
                      operations=operations,
                      fixed_placement=fixed_placement,
                      overlap_mode=None,
                      anchor=None,
                      roi_filter=None,
                      blacklist_filter=None,
                      info=info,
                      genotype=genotype,
                      config_descr=f'vcf_record:{vcf_rec}')
    # end: def import_sv_from_vcf_rec(self, vcf_rec, sim_settings) -> Optional[SV]

    def restore_parent_svs(self, svs: list[SV]) -> list[SV]:
        svs_new: list[SV] = []
        parent2components: dict[str, list[SV]] = defaultdict(list)
        for sv in svs:
            if 'PARENT_SVID' not in sv.info:
                svs_new.append(sv)
            else:
                parent2components[sv.info['PARENT_SVID']].append(sv)

        for parent, components in parent2components.items():
            breakends = []
            for component_sv in components:
                assert component_sv.fixed_placement is not None
                for breakend in component_sv.breakends:
                    breakends.append((component_sv.fixed_placement[breakend], component_sv.sv_id, breakend))
            breakends.sort()
            breakend_map: dict[str, dict[Breakend, Breakend]] = defaultdict(dict)
            parent_placement: list[Locus] = []
            for locus, component_svid, component_breakend in breakends:
                if not parent_placement or parent_placement[-1] != locus:
                    parent_placement.append(locus)
                breakend_map[component_svid][component_breakend] = len(parent_placement)-1
            parent_operations: list[Operation] = []
            for component_sv in components:
                def get_parent_breakend(component_breakend: Breakend) -> Breakend:
                    return breakend_map[component_sv.sv_id][component_breakend]
                for operation in component_sv.operations:
                    parent_operation = Operation(transform=operation.transform,
                                                 novel_insertion_seq=operation.novel_insertion_seq)
                    if operation.source_breakend_region is not None:
                        parent_operation.source_breakend_region = BreakendRegion(
                            get_parent_breakend(operation.source_breakend_region.start_breakend),
                            get_parent_breakend(operation.source_breakend_region.end_breakend))
                    if operation.target_insertion_breakend is not None:
                        assert operation.target_insertion_order is not None
                        parent_operation.target_insertion_breakend = get_parent_breakend(operation.target_insertion_breakend)
                        parent_operation.target_insertion_order = operation.target_insertion_order[1:]
                    parent_operations.append(parent_operation)
                if component_sv.is_interchromosomal:
                    is_interchromosomal = True

            breakend_interval_lengths = [
                locus2.pos - locus1.pos if locus1.chrom == locus2.chrom else None
                for locus1, locus2 in pairwise(parent_placement)]

            info = dict(components[0].info)
            info['SVTYPE'] = info.pop('PARENT_SVTYPE')
            
            svs_new.append(BaseSV(
                sv_id=self.make_sv_id(),
                breakend_interval_lengths=breakend_interval_lengths,
                breakend_interval_min_lengths=[None] * len(breakend_interval_lengths),
                is_interchromosomal=any(component_sv.is_interchromosomal for component_sv in components),
                operations=parent_operations,
                fixed_placement=parent_placement,
                overlap_mode=None,
                anchor=None,
                roi_filter=None,
                blacklist_filter=None,
                info=info,
                genotype=components[0].genotype,
                config_descr=','.join(component_sv.config_descr for component_sv in components)))
        return svs_new
    # end: def restore_parent_svs(self, svs: list[SV]) -> list[SV]

# end: class ImportedVariantSetMaker(VariantSetMaker)

#############################################

VARIANT_SET_MAKER_CLASSES: list[Type[VariantSetMaker]] = [
    FromGrammarVariantSetMaker,
    TandemRepeatVariantSetMaker,
    ChromosomalTranslocationVariantSetMaker,
    ImportedVariantSetMaker
]

def make_variant_set_from_config(vset_config, sim_settings) -> list[SV]:  # type: ignore
    with error_context(vset_config):
        for variant_set_maker_class in VARIANT_SET_MAKER_CLASSES:
            if variant_set_maker_class.can_make_from(vset_config):
                return variant_set_maker_class(vset_config, sim_settings).make_variant_set()
        chk(False, f"Do not know how to make variant set from given config")

VCF_HEADER_INFOS=[dict(id='END', number=1, type='Integer', description="End position of the variant described in this record"),
                  dict(id='INSSEQ', number=1, type='String', description="Novel insertion sequence"),
                  dict(id='INSORD', number=1, type='Integer',
                       description="Insertion order for insertions at given position"),
                  dict(id='VSET', number=1, type='Integer',
                       description="Variant set number (numbered from 0) from which this variant was created"),
                  dict(id='SVTYPE', number=1, type='String',
                       description="Type of structural variant"),
                  dict(id='SVLEN', number=1, type='Integer',
                       description="Length of structural variant"),
                  dict(id='NCOPIES', number=1, type='Integer',
                       description="Number of sequence copies to insert at target"),
                  dict(id='TARGET', number=1, type='Integer',
                       description="Target location for a dispersed duplication or translocation"),
                  dict(id='TARGET_CHROM', number=1, type='String',
                       description="Target chromosome for a dispersed duplication or translocation"),
                  dict(id='OVLP', number=1, type='String',
                       description="Type of overlap region on which the SV or its parent was placed"),
                  dict(id='PARENT_SVID', number=1, type='String',
                       description="ID of parent SV of which this record is one part"),
                  dict(id='PARENT_SVTYPE', number=1, type='String',
                       description="type of parent SV of which this record is one part")]


def get_vcf_header_infos() -> list[dict]:
    vcf_header_infos = copy.deepcopy(VCF_HEADER_INFOS)
    for variant_set_maker_class in VARIANT_SET_MAKER_CLASSES:
        vcf_header_infos.extend(variant_set_maker_class.get_vcf_header_infos())
    return vcf_header_infos

#############################################
