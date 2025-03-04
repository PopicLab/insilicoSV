from abc import ABC, abstractmethod
from collections import defaultdict
from contextlib import closing
import copy
import math
import logging
import random
from typing_extensions import ClassVar, Type, override
import re

from pysam import FastaFile, VariantFile

from insilicosv import utils
from insilicosv.utils import RegionFilter, OverlapMode, Locus, error_context, chk, TandemRepeatRegionFilter
from insilicosv.sv_defs import (Transform, TransformType, BreakendRegion, Operation, SV, VariantType, BaseSV,
                                Syntax, Symbol, SV_KEY, TandemRepeatExpansionContractionSV)

logger = logging.getLogger(__name__)

class VariantSet(ABC):

    @classmethod
    @abstractmethod
    def can_make_from(cls, vset_config):
        raise NotImplementedError()

    @abstractmethod
    def make_variant_set(self):
        raise NotImplementedError()

    @classmethod
    def get_vcf_header_infos(cls):
        return []

    next_sv_id: ClassVar[int] = 0

    @staticmethod
    def make_sv_id():
        sv_id = f'sv{VariantSet.next_sv_id}'
        VariantSet.next_sv_id += 1
        return sv_id

    def __init__(self, vset_config, config):
        self.vset_config = copy.deepcopy(vset_config)
        self.config = config
        self.overlap_kinds = utils.as_list(self.vset_config.get('overlap_region_type', 'all'))
        self.overlap_ranges = []
        self.header = []
        self.overlap_mode = None
        if 'overlap_mode' in self.vset_config:
            chk(isinstance(self.vset_config['overlap_mode'], str) or
                (isinstance(self.vset_config['overlap_mode'], list) and all(isinstance(ovlp, str) for ovlp in self.vset_config['overlap_mode'])),
                'overlap_mode not a string or list of strings')
            try:
                self.overlap_mode = OverlapMode(self.vset_config['overlap_mode'])
            except ValueError:
                chk(False, 'Invalid overlap_mode')
        self.novel_insertion_seqs = None

        self.vset_config['config_descr'] = ', '.join("%s: %s" % item for item in self.vset_config.items())

        chk(isinstance(self.vset_config, dict), f'Each variant set must be a dict')

        self.vset_config['overlap_region_type'] = (tuple(utils.as_list(self.vset_config['overlap_region_type']))
                                                   if 'overlap_region_type' in self.vset_config else ('all',))
        if self.overlap_mode is not None and 'overlap_region_type' not in self.vset_config:
            logger.warning(f'overlap_region_type is not specified for {self} all provided ROIs will be considered')
        self.overlap_ranges = tuple(
            self.vset_config['overlap_region_length_range'] if 'overlap_region_length_range' in self.vset_config else
            [None, None])


    def get_sampled_int_value(self, value, locals_dict=
    None):
        if isinstance(value, (int, float)):
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

    def grammar_to_variant_set(self, lhs_strs, rhs_strs, symbol_lengths, symbol_min_lengths, num_letters, novel_insertion_seqs,
                               n_copies_list, divergence_prob_list, replacement_seq=None):
        #
        # Parse the LHS strings into Symbols, and RHS strings into RHSItems.
        # Create unique Symbols for dispersions.
        #
        # If the overlap anchor was specified, determine the start and end
        # breakends of the anchor.
        #
        overlap_anchor_bounds = [None, None]
        operations = []
        # Retains the symbols representing novel insertions and their corresponding novel sequence.
        novel_insertions: dict[Symbol, str] = {}
        lhs: list[Symbol] = []
        lhs_dispersion: list[int] = []
        breakend_interval_lengths = []
        breakend_interval_min_lengths = []

        # Parse the left hand side of the  grammar, recover the symbols and the potential anchor.
        # The lhs is used to defined the source breakend regions. Because there is the same number of dispersion in the
        # right and left hand side, no new breakend can appear in the rhs (target insertion breakend will be on source breakends from the lhs).
        for lhs_str in lhs_strs:
            if lhs_str == Syntax.ANCHOR_START:
                overlap_anchor_bounds[0] = len(lhs)
                continue

            if lhs_str == Syntax.ANCHOR_END:
                overlap_anchor_bounds[1] = len(lhs)
                continue

            if lhs_str == Syntax.DISPERSION:
                chk(len(symbol_lengths) >= num_letters + len(lhs_dispersion) + 1, f'Invalid dispersion ranges for {lhs_str}')
                breakend_interval_lengths.append(symbol_lengths[num_letters + len(lhs_dispersion)])
                breakend_interval_min_lengths.append(symbol_min_lengths[num_letters + len(lhs_dispersion)])
                lhs_dispersion.append(len(lhs))
            else:
                chk((len(lhs_str) == 1 and lhs_str[0].isupper()),f'Invalid LHS symbol for {lhs_str}')
                chk(len(symbol_lengths) >= len(lhs)-len(lhs_dispersion), f'Missing length ranges {lhs_str}')
                breakend_interval_lengths.append(symbol_lengths[len(lhs)-len(lhs_dispersion)])
                breakend_interval_min_lengths.append(symbol_min_lengths[len(lhs)-len(lhs_dispersion)])
            lhs.append(Symbol(lhs_str))
        chk((overlap_anchor_bounds[0] is None) or (overlap_anchor_bounds[1] is not None), f'Anchor opening without closing: {lhs_strs}')
        overlap_anchor = None
        if overlap_anchor_bounds[0] is not None and overlap_anchor_bounds[1] is not None:
            overlap_anchor = BreakendRegion(overlap_anchor_bounds[0], overlap_anchor_bounds[1])

        n_dispersions_rhs = 0
        current_breakend = 0
        current_insertion_order = 0
        n_multiple_copies = 0
        n_divergence_prob = 0
        # Keep track of the letters that have been deleted or moved for which DEL operations have to be added.
        delete_letters = {letter: True for letter in lhs if letter.name != Syntax.DISPERSION}
        # Parse the rhs to define the operations.
        for rhs_str in rhs_strs:
            chk(rhs_str not in [Syntax.ANCHOR_START, Syntax.ANCHOR_END], f'The anchor can only be specified in the '
                                                                         f'source. If for instance, you want to constrain '
                                                                         f'the insertion target of a dDUP, '
                                                                         f'write A_(). {rhs_strs}')
            if rhs_str == Syntax.DISPERSION:
                n_dispersions_rhs += 1
                chk(n_dispersions_rhs - 1 < len(lhs_dispersion),
                    f'The number of dispersions needs to be the same in the left and right hand sides. {lhs_strs} -> {rhs_strs}')
                # The dispersions are in place and used to determine what moved.
                # When finding a dispersion in the rhs we move the current breakend
                # after the right breakend of the dispersion.
                current_breakend = lhs_dispersion[n_dispersions_rhs - 1] + 1
                continue
            symbol = Symbol(rhs_str[0].upper())
            # The operation is not in place if the symbol does not appear in the lhs or there is a "+".
            is_in_place = symbol in lhs
            if is_in_place:
                source_position = lhs.index(symbol)
                # If there is a different number of dispersions before the symbol in the lhs and rhs then the symbol
                # has been moved.
                is_in_place *= n_dispersions_rhs == len([disp for disp in lhs_dispersion if disp < source_position])
                # If it is still in place we check if the current_breakend is higher or lower to the symbol's
                # source breakend. If it is lower we arbitrarily say it is in place. For contiguous letters it is equivalent.
                is_in_place *= current_breakend <= source_position
            if is_in_place:
                current_breakend = source_position + 1
                # This letter has an inplace version, it should not be deleted.
                delete_letters[symbol] = False

            transform_type = TransformType.IDENTITY
            if rhs_str[0].islower():
                transform_type = TransformType.INV

            # Determine the number of copies of a symbol are needed (when "+" appears in the rhs)
            n_copies = 1
            if Syntax.MULTIPLE_COPIES in rhs_str:
                chk(n_multiple_copies < len(n_copies_list), f'A number of copies must be provided '
                                                                              f'for each + symbol used.')
                n_copies = self.get_sampled_int_value(n_copies_list[n_multiple_copies])
                chk(n_copies >= 1, f'The number of copies must be strictly positive, got {n_copies} for '
                                   f'the {n_multiple_copies + 1} \'+\' symbol {symbol}.')
                n_multiple_copies += 1

            divergence_prob = 0
            if Syntax.DIVERGENCE in rhs_str and replacement_seq is None:
                # replacement_seq is used to provide the SNP when loading from vcf
                chk(n_divergence_prob < len(divergence_prob_list), f'A number of copies must be provided '
                                                            f'for each + symbol used.')
                divergence_prob = self.get_sampled_int_value(divergence_prob_list[n_divergence_prob])
                chk(0 < divergence_prob <= 1, f'The divergence probability must be between 0 (excluded) and 1 (included), got {divergence_prob} for '
                                   f'the {n_divergence_prob + 1} \'*\' symbol {symbol}.')
                n_divergence_prob += 1


            transform = Transform(
                transform_type=transform_type,
                is_in_place=is_in_place,
                divergence_prob=divergence_prob,
                n_copies=n_copies,
                replacement_seq=replacement_seq
            )
            # We do not add operations for inplace identity transformations without divergence (do not affect the sequence).
            if (transform_type != TransformType.IDENTITY or not is_in_place or divergence_prob > 0 or replacement_seq is not None
                    or n_copies > 1):
                operation = Operation(transform=transform, op_info={'SOURCE_LETTER': symbol.name})
                if is_in_place and not n_copies > 1:
                    operation.source_breakend_region = BreakendRegion(current_breakend - 1, current_breakend)
                else:
                    # The insertion order helps to determine how to order multiple events inserted at the same breakend.
                    operation.target_insertion_order = (current_insertion_order,)
                    current_insertion_order += 1
                    operation.target_insertion_breakend = current_breakend
                    if symbol in lhs:
                        operation.source_breakend_region = BreakendRegion(source_position, source_position + 1)
                    else:
                        if symbol not in novel_insertions:
                            # The position of the length of a new letter is after all the letters in the lhs (non dispersion)
                            # and in order of apparition in the rhs
                            chk(len(symbol_lengths) >= len(lhs)-len(lhs_dispersion)+len(novel_insertions),
                                f'missing length_ranges for {symbol} in {lhs_strs}')
                            novel_insertion = None
                            if novel_insertion_seqs is not None:
                                chk(len(novel_insertion_seqs) > len(novel_insertions), f'If a novel_insertions field is provided, it must '
                                                                                        f'contain the sequences for each novel insertion or None.')
                                if novel_insertion_seqs[len(novel_insertions)] is not None:
                                    novel_insertion = novel_insertion_seqs[len(novel_insertions)]
                            if novel_insertion is None:
                                novel_insertion = utils.generate_seq(length=symbol_lengths[len(lhs)-len(lhs_dispersion)+len(novel_insertions)])
                            novel_insertions[symbol] = novel_insertion
                        operation.novel_insertion_seq = novel_insertions[symbol]
                operations.append(operation)
        for letter, delete in delete_letters.items():
            if delete:
                source_index = lhs.index(letter)
                # The symbol was deleted, we add a DEL operation.
                transform = Transform(transform_type=TransformType.DEL, is_in_place=True)
                operation = Operation(transform=transform, op_info={'SOURCE_LETTER': letter.name})
                operation.source_breakend_region = BreakendRegion(source_index, source_index + 1)
                operations.append(operation)
        chk(len(lhs_dispersion) == n_dispersions_rhs, f'Dispersions cannot be altered, the same dispersions must '
                                                      f'be present on the left and right sides. lhs: {len(lhs_dispersion)} dispersions,'
                                                      f'rhs: {n_dispersions_rhs} dispersions.')
        # If overlap_anchor is not specified, infer "full SV" as the anchor if no dispersion
        if (overlap_anchor is None) and (self.overlap_mode is not None) and (not lhs_dispersion):
            overlap_anchor = BreakendRegion(0, len(lhs)-1)
        return operations, overlap_anchor, lhs_dispersion, breakend_interval_lengths, breakend_interval_min_lengths

    # end: def grammar_to_variant_set


class SimulatedVariantSet(VariantSet):

    def __init__(self, vset_config, config):
        super().__init__(vset_config, config)

        chk(isinstance(vset_config.get('type'), str), f'Missing or bad variant type')
        chk(isinstance(vset_config.get('number'), int) and vset_config['number'] >= 0,
            f'Missing or bad "number" of variants to generate')

        try:
            self.vset_config['type'] = VariantType(vset_config['type'])
        except ValueError:
            chk(False, 'Invalid variant type')

        self.preprocess_config()

    def preprocess_config(self):
        chk('divergence_prob' not in self.vset_config or isinstance(self.vset_config['divergence_prob'], list),
            'divergence_prob must be a list of floats in ]0, 1] or ranges')

        if (self.vset_config['type'] == VariantType.DUP) and ('n_copies' not in self.vset_config):
            self.vset_config['n_copies'] = [1]
        if 'n_copies' in self.vset_config:
            chk(isinstance(self.vset_config['n_copies'], list), 'The number of copies must a list of integers or ranges.')
        if self.vset_config['type'] == "mCNV":
            chk(('n_copies' in self.vset_config) and (self.vset_config['n_copies'][0] > 1),
                f'n_copies has to be provided and be above 1 for a mCNV.')
    # end: def preprocess_config(self)
        
    @property
    def sv_type(self):
        return self.vset_config['type']

    def get_roi_filter(self):
        if self.overlap_mode is None:
            return None
        return RegionFilter(region_kinds=self.overlap_kinds, region_length_range=self.overlap_ranges)

    def get_blacklist_filter(self):
        vset_config = self.vset_config
        if 'blacklist_region_type' not in vset_config:
            return None
        return RegionFilter(region_kinds=tuple(utils.as_list(self.vset_config['blacklist_region_type'])))


    def pick_genotype(self):
        if self.config.get('homozygous_only', False) or random.randint(0, 1):
            return True, True
        else:
            return random.choice([(True, False), (False, True)])

    @override
    def make_variant_set(self):
        return [self.simulate_sv() for _ in range(self.vset_config['number'])]

    @abstractmethod
    def simulate_sv(self):
        raise NotImplementedError()

#############################################
# Constructing SVs from grammar definitions #
#############################################

class FromGrammarVariantSet(SimulatedVariantSet):
    """Constructs an SV from a grammar definition"""

    @override
    @classmethod
    def can_make_from(cls, vset_config):
        return vset_config.get("type") in [k.value for k in SV_KEY.keys()] + ['Custom']

    @override
    def preprocess_config(self):
        super().preprocess_config()

        for vset_config_key in self.vset_config:
            chk(vset_config_key in (
                'type', 'number',
                'source', 'target',
                'length_ranges', 'dispersion_ranges',
                'overlap_region_type', 'overlap_region_length_range',
                'overlap_mode',
                'blacklist_region_type',
                'divergence_prob',
                'n_copies',
                'novel_insertions',
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
            chk('divergence_prob' not in vset_cfg or vset_cfg['divergence_prob'] == [1],
                'divergence prob for SNP can only be 1')
            vset_cfg['length_ranges'] = [[1, 1]]
            vset_cfg['divergence_prob'] = [1.0]

        if 'novel_insertions' in self.vset_config:
            try:
                with open(self.vset_config['novel_insertions'], 'r') as sequences:
                    self.novel_insertion_seqs = [line.rstrip() for line in sequences]
                    chk(all(bool(re.match('^[TCGA]+$', line)) for line in self.novel_insertion_seqs), f'The file novel_insertions'
                                                                                  f'can only contains invalid characters. It must be a list of sequences.')
            except:
                chk(False, f'novel_insertion must be a readable file containing a sequence per line.')
        for src_trg in ('source', 'target'):
            chk(isinstance(self.vset_config.get(src_trg), (type(None), list, str, tuple)),
                f'{src_trg} must be a string or list of strings')
            chk(not isinstance(self.vset_config.get(src_trg), list) or
                utils.is_list_of(str, self.vset_config.get(src_trg)),
                f'elements of {src_trg} must be strings')

        if self.vset_config['type']  == VariantType.Custom:
            chk('source' in self.vset_config and 'target' in self.vset_config,
                'Please provide "source" and "target" grammar for custom SV')
        else:
            lhs_strs, rhs_strs = SV_KEY[VariantType(self.vset_config['type'] )]
            if 'source' in self.vset_config:
                lhs_strs = tuple(self.vset_config['source'])
                chk(tuple(s for s in lhs_strs if s not in (Syntax.ANCHOR_START, Syntax.ANCHOR_END)) ==
                    SV_KEY[VariantType(self.vset_config['type'] )][0],
                    f'Source spec does not match built-in spec')
            else:
                self.vset_config['source'] = lhs_strs
            chk('target' not in self.vset_config, f'Please do not provide a target for predefined types.'
                                                  f'If needed, anchors must be specified in the source')
            self.vset_config['target'] = rhs_strs

        rhs_strs_list: list[str] = []
        for c in self.vset_config['target']:
            if c in (Syntax.DIVERGENCE, Syntax.MULTIPLE_COPIES):
                chk(rhs_strs_list and rhs_strs_list[-1] and rhs_strs_list[-1][0].isalpha(),
                    f'{c} must modify a symbol')
                rhs_strs_list[-1] += c
            else:
                rhs_strs_list.append(c)

        self.vset_config['target'] = rhs_strs_list

        # Anchor the whole SV if overlap mode is specified, no anchor is provided and there is no dispersion.
        if ((self.overlap_mode is not None) and (Syntax.DISPERSION not in self.vset_config['source']) and
                (Syntax.ANCHOR_START not in self.vset_config['source'])):
            self.vset_config['source'] = tuple([Syntax.ANCHOR_START, *self.vset_config['source'], Syntax.ANCHOR_END])

    @property
    def is_interchromosomal(self):
        return self.vset_config.get('interchromosomal', False)

    def symmetrize(self):
        lhs_strs = self.vset_config['source']
        rhs_strs = self.vset_config['target']
        # Enforce the symmetry of the predefined SVs with duplications or dispersions.
        if (("DUP" in self.vset_config['type'].name or "TRA" in self.vset_config['type'].name or
                 "iDEL" in self.vset_config['type'].name)
                 and (not self.is_interchromosomal)
                 and random.randint(0, 1)):
            def flip_anchor(val: str) -> str:
                if val == Syntax.ANCHOR_START:
                    return Syntax.ANCHOR_END
                if val == Syntax.ANCHOR_END:
                    return Syntax.ANCHOR_START
                return val
            lhs_strs = tuple(map(flip_anchor, lhs_strs))[::-1]
            rhs_strs = tuple(map(flip_anchor, rhs_strs))[::-1]
            if "length_ranges" in self.vset_config:
                self.vset_config['length_ranges'] = self.vset_config['length_ranges'][::-1]
            if "dispersion_ranges" in self.vset_config:
                self.vset_config['dispersion_ranges'] = self.vset_config['dispersion_ranges'][::-1]
        chk(all(len(rhs)<3 for rhs in rhs_strs), f'The operators + and * cannot be used at the same time {rhs_strs}.')
        return lhs_strs, rhs_strs

    # Randomly pick distances withing the ranges define for each symbol
    @staticmethod
    def pick_symbol_lengths(length_ranges, dispersion_ranges, letter_indexes):
        ranges = length_ranges + dispersion_ranges
        remaining_symbols = [i for i in range(len(ranges))]
        chk(all(isinstance(length_range, str) or (isinstance(length_range, list) and (len(length_range) == 2))
                for length_range in ranges),f'length_ranges must be a list of [min, max] pairs')
        # Keep track of the dependencies between the symbols length ranges.
        # A dependency is the index of the letter the range is depending on (possibly a different letter for the min and max bounds)
        # and the offset to those letter lengths.
        dependencies = {}
        symbol_lengths = {}
        symbol_min_lengths = {}
        while remaining_symbols:
            idx = remaining_symbols.pop(0)
            min_range, max_range = ranges[idx]
            if(isinstance(min_range, int) or min_range is None) and (isinstance(max_range, int) or max_range is None):
                # Both bounds are independent to other letter lengths.
                def assign_length(min_range, max_range, is_dispersion):
                    if max_range is None:
                        length = None
                        if min_range is not None:
                            chk(is_dispersion, f'Only dispersions can have min but not max length')
                        min_length = min_range
                    else:
                        chk(min_range is not None, f'max_length given but not min_length')
                        chk(min_range <= max_range, f'max bound less than min bound')
                        chk(min_range >= 0, 'min length cannot be negative')
                        length = random.randint(min_range, max_range)
                        min_length = None
                    return length, min_length

                length, min_length = assign_length(min_range, max_range, idx >= len(length_ranges))
                symbol_lengths[idx] = length
                symbol_min_lengths[idx] = min_length
            else:
                for pos, bound in enumerate([min_range, max_range]):
                    if isinstance(bound, str):
                        # The bound is dependent of another letter
                        no_space_bound = bound.strip()
                        format_bound = [int(char) if char.isdigit() else char for char in no_space_bound]
                        letters_bound= [(idx, char) for idx, char in enumerate(format_bound) if isinstance(char, str) and char.isalpha()]
                        chk(all(letter in letter_indexes for _, letter in letters_bound), f'The length of a symbol depends on {letters_bound}'
                                                                                    f'one of which is not define in neither the source nor target.')
                        computed_letter = 0
                        for idx_in_dependency, letter in letters_bound:
                            # Gets the index of the letter in the list of length ranges.
                            index = letter_indexes[letter]
                            if index in symbol_lengths:
                                chk(symbol_lengths[index] is not None, f'A symbol length cannot depend on an unbounded symbol.')
                                ope = ''
                                if (idx_in_dependency > 0) and (isinstance(format_bound[idx_in_dependency-1], int)
                                                                or format_bound[idx_in_dependency-1] in letter_indexes):
                                    # The previous character was an integer, so a multiplication symbol was omitted.
                                    ope = '*'
                                format_bound[idx_in_dependency] = ope + str(symbol_lengths[index])
                                computed_letter += 1
                            else:
                                if not idx in dependencies:
                                    dependencies[idx] = []
                                dependencies[idx].append(index)
                                if index in dependencies:
                                    chk(not idx in dependencies[index], f'There is a cyclic dependency in the length definitions {letter}.')
                                    dependencies[idx] += dependencies[index]
                                dependencies[idx] = list(set(dependencies[idx]))
                        if len(letters_bound) == computed_letter:
                            # All the dependencies are satisfied, we evaluate the operation
                            try:
                                # We try to update the bound if the operation is valid
                                eval_bound = eval("".join([str(char) for char in format_bound]))
                                if pos == 0:
                                    eval_bound = math.ceil(eval_bound)
                                else:
                                    eval_bound = math.floor(eval_bound)
                                ranges[idx][pos] = eval_bound
                            except:
                                chk(False, f'The length of a symbol depends on an invalid operation {idx}.')
                            chk(ranges[idx][pos] >= 0, f'The operation defining the {idx}th symbol length gives a negative bound.')
                remaining_symbols.append(idx)
        lengths = sorted(symbol_lengths.items(), key=lambda x: x[0])
        min_lengths = sorted(symbol_min_lengths.items(), key=lambda x: x[0])
        return [l[1] for l in lengths], [m_l[1] for m_l in min_lengths]

    @override
    def simulate_sv(self) -> SV:
        lhs_strs, rhs_strs = self.symmetrize()
        length_ranges = self.vset_config['length_ranges'] if 'length_ranges' in self.vset_config else []
        dispersion_ranges = self.vset_config['dispersion_ranges'] if 'dispersion_ranges' in self.vset_config else []
        num_dispersion = len([letter for letter in lhs_strs if letter == Syntax.DISPERSION])
        chk(num_dispersion == len(dispersion_ranges),
            f'Missing dispersion length ranges, {len(dispersion_ranges)} provided, {num_dispersion} expected')
        letters = [letter for letter in lhs_strs if letter not in [Syntax.ANCHOR_END, Syntax.ANCHOR_START, Syntax.DISPERSION]]
        chk(len(letters) == len(set(letters)), f'Duplicate LHS symbol {letters}')

        # Add novel insertion letters only appearing in the rhs.
        for letter in rhs_strs:
            if letter[0].upper() not in letters + [Syntax.ANCHOR_END, Syntax.ANCHOR_START, Syntax.DISPERSION]:
                letters.append(letter[0].upper())
        chk(len(length_ranges) == len(letters), f'Missing length ranges, expected {len(letters)} provided {len(length_ranges)}')
        letter_indexes = {letter: index for index, letter in enumerate(letters)}

        # Compute the lengths of the different symbols and dispersions from the length ranges.
        symbol_lengths, symbol_min_lengths = self.pick_symbol_lengths(length_ranges, dispersion_ranges, letter_indexes)

        novel_insertion_seqs = self.novel_insertion_seqs
        n_copies_list = self.vset_config.get('n_copies', [])
        divergence_prob_list = self.vset_config.get('divergence_prob', [])

        # Build the different operations and anchor, determine the breakends and the distance between them.
        (operations, anchor, dispersions, breakend_interval_lengths, 
         breakend_interval_min_lengths) = self.grammar_to_variant_set(lhs_strs, rhs_strs, symbol_lengths, symbol_min_lengths, len(letters),
                                                                      novel_insertion_seqs, n_copies_list, divergence_prob_list)
        #
        # construct the SV object
        #
        info = self.construct_info()
        roi_filter = self.get_roi_filter()
        return BaseSV(sv_id=self.make_sv_id(),
                      breakend_interval_lengths=breakend_interval_lengths,
                      breakend_interval_min_lengths=breakend_interval_min_lengths,
                      is_interchromosomal=self.is_interchromosomal,
                      operations=operations,
                      anchor=anchor,
                      dispersions=dispersions,
                      overlap_mode=self.overlap_mode,
                      roi_filter=roi_filter,
                      blacklist_filter=self.get_blacklist_filter(),
                      fixed_placement=None,
                      info=info,
                      genotype=self.pick_genotype(),
                      config_descr=self.vset_config['config_descr'])

    def construct_info(self):
        sv_type_str = self.vset_config['type'] .value
        source_str = ''.join(self.vset_config['source'])
        target_str = ''.join(self.vset_config['target'])
        grammar = f'{source_str}>{target_str}'
        return {'SVTYPE': sv_type_str, 'GRAMMAR': grammar}

# end: class FromGrammarVariantSetMaker

class TandemRepeatVariantSet(SimulatedVariantSet):

    @override
    @classmethod
    def can_make_from(cls, vset_config):
        return vset_config.get("type") in ('trEXP', 'trCON')

    def __init__(self, vset_config, config):
        super().__init__(vset_config, config)
        chk('repeat_count_change_range' in self.vset_config,
            'Please specify repeat_count_change_range')
        chk(utils.is_valid_closed_int_range(self.vset_config['repeat_count_change_range']),
            'Invalid repeat_count_change_range')

    @override
    def preprocess_config(self):
        for vset_config_key in self.vset_config:
            chk(vset_config_key in (
                'type', 'number',
                'repeat_count_change_range',
                'overlap_region_type', 'overlap_region_length_range',
                'overlap_mode',
                'blacklist_region_type',
                'config_descr',
            ), f'invalid SV config key {vset_config_key}')

    @override
    def simulate_sv(self):

        repeat_count_change = random.randint(*self.vset_config['repeat_count_change_range'])
        info = dict(SVTYPE=self.vset_config['type'].value, TR_CHANGE=repeat_count_change)
        overlap_region_type = (tuple(utils.as_list(self.vset_config['overlap_region_type']))
                               if 'overlap_region_type' in self.vset_config else 'all')
        anchor = BreakendRegion(0, 1)
        self.overlap_mode = OverlapMode.EXACT

        chk(self.overlap_mode in [None, OverlapMode.EXACT], 'overlap_mode must be "exact" for Tandem Repeat')
        if self.vset_config['type'] == VariantType.trEXP:
            breakend_interval_lengths = [None]
            operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY,
                                                        is_in_place=False,
                                                        n_copies=repeat_count_change),
                                    source_breakend_region=BreakendRegion(start_breakend=0, end_breakend=1),
                                    target_insertion_breakend=0,
                                    target_insertion_order=(0,),
                                    op_info={'SOURCE_LETTER': 'A'})]
            # We only need the repeat motif to be present once.
            roi_filter = TandemRepeatRegionFilter(min_num_repeats=1,
                                                  region_kinds=overlap_region_type)
            return TandemRepeatExpansionContractionSV(
                sv_id=self.make_sv_id(),
                breakend_interval_lengths=breakend_interval_lengths,
                breakend_interval_min_lengths=[None] * len(breakend_interval_lengths),
                is_interchromosomal=False,
                operations=operations,
                anchor=anchor,
                overlap_mode=self.overlap_mode,
                roi_filter=roi_filter,
                blacklist_filter=self.get_blacklist_filter(),
                fixed_placement=None,
                info=info,
                genotype=self.pick_genotype(),
                config_descr=self.vset_config['config_descr'],
                num_repeats_in_placement=1,
                dispersions=[])
        elif self.vset_config['type'] == VariantType.trCON:
            breakend_interval_lengths = [None]
            operations = [Operation(transform=Transform(transform_type=TransformType.DEL,
                                                        is_in_place=True, n_copies=1),
                                    source_breakend_region=BreakendRegion(0, 1),
                                    op_info={'SOURCE_LETTER': 'A'})]

            # ensure there are enough existing repeats to delete.
            roi_filter = TandemRepeatRegionFilter(min_num_repeats=repeat_count_change,
                                                  region_kinds=overlap_region_type)

            return TandemRepeatExpansionContractionSV(
                sv_id=self.make_sv_id(),
                breakend_interval_lengths=breakend_interval_lengths,
                breakend_interval_min_lengths=[None],
                is_interchromosomal=False,
                operations=operations,
                anchor=anchor,
                overlap_mode=self.overlap_mode,
                roi_filter=roi_filter,
                blacklist_filter=self.get_blacklist_filter(),
                fixed_placement=None,
                info=info,
                genotype=self.pick_genotype(),
                config_descr=self.vset_config['config_descr'],
                num_repeats_in_placement=repeat_count_change,
                dispersions=[])
        else:
            assert False

    @override
    @classmethod
    def get_vcf_header_infos(cls):
        return [dict(id='TR_INS_REP_SEQ', number=1, type='String',
                            description="for tandem repeat insertions, the sequence of a single repeat"),
                dict(id='TR_CHANGE', number=1, type='Integer',
                     description="for tandem repeat insertions/expansions/deletions, number of "
                     "repeats inserted or removed")]

#############################################

###########################
# Importing SVs from VCFs #
###########################

class ImportedVariantSet(VariantSet):

    @override
    @classmethod
    def can_make_from(cls, vset_config):
        return 'import' in vset_config

    def __init__(self, vset_config, config):
        super().__init__(vset_config, config)
        chk(utils.is_readable_file(vset_config['import']),
            '{path} vcf must name a readable file'.format(path=vset_config['import']))
        chk(set(vset_config.keys()) <= {'import'}, 'invalid config key')

        with FastaFile(config['reference']) as reference:
            self.chrom_lengths = {chrom: chrom_length
                                  for chrom, chrom_length in zip(reference.references,
                                                                 reference.lengths)}

    @override
    def make_variant_set(self):
        recs = defaultdict(list)
        num_simple_sv = 0
        with closing(VariantFile(self.vset_config['import'])) as vcf:
            self.header = vcf.header
            for vcf_rec in vcf.fetch():
                with error_context(vcf_rec):
                    if vcf_rec.chrom not in self.chrom_lengths: continue
                    vcf_info = dict(vcf_rec.info)
                    if 'PARENT_SVID' in vcf_info:
                        recs[vcf_info['PARENT_SVID']].append(vcf_rec)
                    else:
                        recs[str(num_simple_sv)].append(vcf_rec)
                        num_simple_sv += 1
        svs = []
        for parent_id, sv_recs in recs.items():
            sv = self.import_sv_from_vcf_recs(sv_recs, parent_id)
            if sv is not None:
                svs.append(sv)
        return svs

    def parse_vcf_rec_info(self, vcf_rec):
        chk(vcf_rec.chrom in self.chrom_lengths, f'unknown contig {vcf_rec.chrom}')
        chk(len(vcf_rec.alleles) == 2, 'Can only import bi-allelic variants from VCF')
        chk(len(vcf_rec.samples) <= 1, 'Can only import VCFs with one sample')
        parsed_info = {}
        if vcf_rec.samples:
            sample = vcf_rec.samples[0]
            if sample['GT'] == (1, 1):
                parsed_info['GENOTYPE'] = (True, True)
            else:
                parsed_info['GENOTYPE'] = random.choice([(True, False), (False, True)])
        else:
            parsed_info['GENOTYPE'] = random.choice([(True, True), (True, False), (False, True)])
        vcf_info = dict(vcf_rec.info)
        if set(''.join(vcf_rec.alleles).upper()) <= set('TCGA'):
            if len(vcf_rec.alleles[0]) == 1 and len(vcf_rec.alleles[1]) == 1:
                # SNP
                vcf_info['SVTYPE'] = 'SNP'
        chk('SVTYPE' in vcf_info, 'Need an SVTYPE to import from vcf records')
        rec_type_str = vcf_info['SVTYPE']
        can_import_types = [v.value for v in VariantType] + ['IDENTITY']
        chk(rec_type_str in {variant_type for variant_type in can_import_types},
            f'Currently only the following VCF types are supported: {can_import_types} but {rec_type_str} was provided')
        if rec_type_str != 'IDENTITY':
            parsed_info['SVTYPE'] = VariantType(rec_type_str)
        else:
            parsed_info['SVTYPE'] = 'IDENTITY'

        chk(vcf_rec.start is not None, f'The Start position must not be None')
        rec_start = Locus(chrom=vcf_rec.chrom, pos=vcf_rec.start)
        parsed_info['START'] = rec_start
        rec_end = Locus(chrom=vcf_rec.chrom, pos=vcf_rec.stop)

        if 'SVLEN' in vcf_info:
            rec_len = vcf_info['SVLEN']
        else:
            rec_len = rec_end.pos - rec_start.pos
        parsed_info['END'] = rec_end
        parsed_info['SVLEN'] = rec_len

        rec_target = None
        is_interchromosomal = False
        if 'TARGET' in vcf_info:
            chk(isinstance(vcf_info.get('TARGET'), int), 'TARGET has to be an int representing a position')
            rec_target = Locus(chrom=vcf_info.get('TARGET_CHROM', vcf_rec.chrom), pos=vcf_info['TARGET'] - 1)
            if rec_target.chrom != vcf_rec.chrom:
                is_interchromosomal = True
            elif rec_end is not None:
                # The target has to be outside of the source region
                chk(rec_target.pos >= rec_end.pos or rec_target.pos <= rec_start.pos,"invalid dispersion target")
        parsed_info['TARGET'] = rec_target
        parsed_info['INTERCHROMOSOMAL'] = is_interchromosomal

        novel_insertion_seq = None
        if 'INSSEQ' in vcf_info:
            chk(isinstance(vcf_info['INSSEQ'], str) and all(bp in 'TCGA' for bp in vcf_info['INSSEQ']),
                'INSSEQ has to be a valid sequence, {} was provided'.format(vcf_info['INSSEQ']))
            novel_insertion_seq = [vcf_info['INSSEQ']]
        elif parsed_info['SVTYPE'] == VariantType.INS:
            # It is a novel sequence insertion but the sequence to insert isn't provided, we generate a random sequence
            logger.warning(
                'SV containing a novel sequence insertion but the sequence was not provided, a random sequence is used')
            novel_insertion_seq = [utils.generate_seq(length=rec_len)]
        parsed_info['INSSEQ'] = novel_insertion_seq

        parsed_info['NCOPIES'] = vcf_info.get('NCOPIES', 1)
        parsed_info['IN_PLACE'] = vcf_info.get('IN_PLACE', 'in_place' if rec_target is None else 'paste') == 'in_place'
        parsed_info['INSORD'] = vcf_info.get('INSORD', None)
        parsed_info['DIVERGENCEPROB'] = [0]
        parsed_info['ALT'] = None
        parsed_info['PARENT_SVTYPE'] = vcf_info.get('PARENT_SVTYPE', 'Custom')
        if parsed_info['SVTYPE'] == VariantType.SNP:
            parsed_info['DIVERGENCEPROB'] = [1]
            if vcf_rec.alts is not None:
                parsed_info['ALT'] = vcf_rec.alts[0]
        additional_info = {}
        valid_info = [header for header_list in VCF_HEADER_INFOS for id, header in header_list.items() if id == 'id']
        for info, val in vcf_info.items():
            if info not in parsed_info and info in valid_info and val is not None and val:
                additional_info[info] = val
        return parsed_info, additional_info

    @staticmethod
    def remove_duplicates(ordered_list):
        seen = set()
        seen_add = seen.add
        return [x for x in ordered_list if not (x in seen or seen_add(x))]

    def import_sv_from_vcf_recs(self, sv_recs, parent_id):
        sv_operations = []
        source_regions = []
        targets = []
        genotype = None
        parent_info = {'PARENT_SVID': parent_id}
        positions_per_rec = []
        current_insord = 0

        # Extract the info of all the records of the SV to determine the breakends and placement.
        for vcf_rec in sv_recs:
            parsed_info, additional_info = self.parse_vcf_rec_info(vcf_rec)
            if genotype is None:
                genotype = parsed_info['GENOTYPE']
            else:
                chk(genotype == parsed_info['GENOTYPE'], f'All the records of a same SV must have the same genotype, '
                                                         f'found {vcf_rec}')
            if len(sv_recs) > 1:
                parent_info['SVTYPE'] = parsed_info['PARENT_SVTYPE']
            else:
                parent_info['SVTYPE'] = str(parsed_info['SVTYPE'])
            target_left = False
            source_regions.append([parsed_info['START'], parsed_info['END']])
            placement = [parsed_info['START'], parsed_info['END']]
            # Gather the length information to parse the grammar and get the operations associated to the record
            symbol_lengths = [parsed_info['END'].pos - parsed_info['START'].pos]
            symbol_min_lengths = [None]
            insord = (0, )
            if parsed_info['TARGET'] is not None:
                # We check the relative position of the target compared to start and end
                if  not parsed_info['INTERCHROMOSOMAL'] and parsed_info['TARGET'] <= parsed_info['START']:
                    placement = [parsed_info['TARGET'], parsed_info['START'], parsed_info['END']]
                    target_left = True
                else:
                    placement = [parsed_info['START'], parsed_info['END'], parsed_info['TARGET']]
                    chk(parsed_info['TARGET'] >= parsed_info['END'], f'The position fo the target has to be outside of the source region,'
                                                                     f'{vcf_rec} has a target between the start and end.')
                if not parsed_info['INTERCHROMOSOMAL']:
                    #The target is an intrachromosomal dispersion
                    dispersion_dist = parsed_info['START'].pos - parsed_info['TARGET'].pos if target_left \
                                                        else parsed_info['TARGET'].pos - parsed_info['END'].pos
                    symbol_lengths = symbol_lengths + [dispersion_dist]
                    symbol_min_lengths.append(None)
                else:
                    # There is an interchromosomal dispersion
                    placement = [parsed_info['START'], parsed_info['END'], parsed_info['TARGET']]
                    symbol_lengths.append([None, None])
                    symbol_min_lengths.append(None)
                if parsed_info['TARGET'] in targets:
                    chk(parsed_info['INSORD'] is not None, 'Two sequences inserted at the locus {}, the INSORD field must be provided for these records to prevent any ambiguity.'.format(parsed_info['INSORD']))
                insord = (parsed_info['INSORD'],) if parsed_info['INSORD'] is not None else (current_insord,)
                current_insord = insord[0] + 1
            targets.append(parsed_info['TARGET'])
            positions_per_rec.append(placement)
            is_in_place = parsed_info['IN_PLACE']
            # The SVTYPE of a record must be a predefined SV type or the identity operation.
            if isinstance(parsed_info['SVTYPE'], VariantType) and parsed_info['SVTYPE'] not in [VariantType.DEL, VariantType.INV]:
                lhs_strs, rhs_strs = SV_KEY[parsed_info['SVTYPE']]
                if target_left:
                    # The symmetrical of the SV grammar has been used
                    lhs_strs = lhs_strs[::-1]
                    rhs_strs = rhs_strs[::-1]
                rhs_strs_list = []
                for c in rhs_strs:
                    if c in (Syntax.DIVERGENCE, Syntax.MULTIPLE_COPIES):
                        rhs_strs_list[-1] += c
                    else:
                        rhs_strs_list.append(c)
                # Get the operations record by record
                operations, _, _, _, _ = self.grammar_to_variant_set(lhs_strs, rhs_strs_list, symbol_lengths, symbol_min_lengths, 1,
                                                                     parsed_info['INSSEQ'], [parsed_info['NCOPIES']],
                                                                     parsed_info['DIVERGENCEPROB'], replacement_seq=parsed_info['ALT'])
                for operation in operations:
                    operation.op_info = additional_info
                    operation.target_insertion_order = insord
            else:
                # The record is an atomic operation of a complex SV (DEL, INV or IDENTITY)
                max_breakend = 1
                min_offset = 1
                if not is_in_place:
                    if len(set(placement)) == 3:
                        # Check if the insertion target is the end breakend or not
                        max_breakend = 2
                    elif target_left:
                        # Target and start are the same breakend
                        min_offset = 0
                start_breakend = 0 if not target_left else min_offset
                end_breakend = 1 if not target_left else max_breakend
                target_breakend = None
                if parsed_info['TARGET'] is not None:
                    target_breakend = 0 if target_left else max_breakend
                source_region = BreakendRegion(start_breakend, end_breakend)
                transform =Transform(TransformType[vcf_rec.info['SVTYPE']], is_in_place=is_in_place, n_copies=parsed_info['NCOPIES'],
                                     divergence_prob=parsed_info['DIVERGENCEPROB'][0], replacement_seq=parsed_info['ALT'])
                operations =[Operation(transform, source_breakend_region=source_region, novel_insertion_seq=parsed_info['INSSEQ'],
                                       target_insertion_breakend=target_breakend,
                                       target_insertion_order=insord, op_info=additional_info)]
            sv_operations.append(operations)
        # If we have more than one record we unify the breakends and their positions through the different operations
        if len(sv_recs) > 1:
            # Combine the information from the different records to determine the sets of breakends of the whole SV
            # Remove duplicate positions and insert the target positions in order to find the positions of all the breakends
            placements = self.remove_duplicates([source[region_idx] for source in source_regions for region_idx in [0, 1]] +
                                                [target for target in targets if target is not None])
            # Get the breakends ordered
            placements = sorted(placements, key=lambda x: x.pos)
            # Unify the breakend of the different operations to have a single SV with well-defined breakends.
            for rec_idx, operation_list in enumerate(sv_operations):
                positions = positions_per_rec[rec_idx]
                # Get the indexes in the placement of the positions of the source regions, which corresponds to the ordered breakends in the complex SV
                breakends = [placements.index(pos) for pos in positions]
                for operation in operation_list:
                    # A single record might correspond ot multiple operations, they all share the same source but one might have a target and not the other
                    # for instance nrTRA is a DEL without a target and an Identity with a target
                    ope_target = breakends[operation.target_insertion_breakend] if operation.target_insertion_breakend is not None else None
                    start_breakend = breakends[operation.source_breakend_region.start_breakend]
                    end_breakend = breakends[operation.source_breakend_region.end_breakend]
                    operation.source_breakend_region = BreakendRegion(start_breakend=start_breakend, end_breakend=end_breakend)
                    operation.target_insertion_breakend = ope_target
            sv_operations = [sv_operation for operation_list in sv_operations for sv_operation in operation_list]
        else:
            # There is a single record
            sv_operations = sv_operations[0]
            placements = positions_per_rec[0]
        sv_id = 'Imported_' + str(parent_id)
        return BaseSV(sv_id=sv_id,
                      breakend_interval_lengths=[None]*(len(placements)-1), # Positions are known, the lengths are not needed
                      breakend_interval_min_lengths=[None]*(len(placements)-1),
                      is_interchromosomal=None, #The target chromosome is already known from the fixed_placement
                      operations=sv_operations,
                      fixed_placement=placements,
                      overlap_mode=None,
                      anchor=None,
                      roi_filter=None,
                      blacklist_filter=None,
                      info=parent_info,
                      genotype=genotype,
                      config_descr=f'vcf_record:{vcf_rec}',
                      dispersions=[])

    # end: def import_sv_from_vcf_rec(self, vcf_rec, config) -> Optional[SV]

# end: class ImportedVariantSetMaker(VariantSetMaker)

VARIANT_SET_CLASSES: list[Type[VariantSet]] = [
    FromGrammarVariantSet,
    TandemRepeatVariantSet,
    ImportedVariantSet
]

def make_variant_set_from_config(vset_config, config) -> list[SV]:  # type: ignore
    with error_context(vset_config, config):
        for variant_set_class in VARIANT_SET_CLASSES:
            if variant_set_class.can_make_from(vset_config):
                variant_set = variant_set_class(vset_config, config)
                return variant_set.make_variant_set(), variant_set.overlap_ranges, variant_set.overlap_kinds, variant_set.overlap_mode, variant_set.header
        chk(False, f"The format of the config or the sv_type is not supported {vset_config}")

VCF_HEADER_INFOS=[dict(id='END', number=1, type='Integer', description="End position of the variant described in this record"),
                  dict(id='INSSEQ', number=1, type='String', description="Novel insertion sequence"),
                  dict(id='INSORD', number=1, type='Integer',
                       description="Insertion order for insertions at given position"),
                  dict(id='VSET', number=1, type='Integer',
                       description="Variant set number (numbered from 0) from which this variant was created"),
                  dict(id='SVTYPE', number=1, type='String',
                       description="Type of structural variant"),
                  dict(id='GRAMMAR', number=1, type='String',
                       description="Grammar of the structural variant"),
                  dict(id='SOURCE_LETTER', number=1, type='String',
                       description="Source letter the record altered"),
                  dict(id='IN_PLACE', number=1, type='String',
                       description="If the operation was performed at the source location or pasted a modified source at a target location"),
                  dict(id='SVLEN', number=1, type='Integer',
                       description="Length of structural variant"),
                  dict(id='NCOPIES', number=1, type='Integer',
                       description="Number of sequence copies to insert at target"),
                  dict(id='TARGET', number=1, type='Integer',
                       description="Target location for a dispersed duplication or translocation"),
                  dict(id='TARGET_CHROM', number=1, type='String',
                       description="Target chromosome for a dispersed duplication or translocation"),
                  dict(id='OVLP', number=1, type='String',
                       description="Type of ROI on which the SV component was placed"),
                    dict(id='OVLP_TARGET', number=1, type='String',
                       description="Type of ROI on which the insertion target of an SV component was placed"),
                  dict(id='PARENT_SVID', number=1, type='String',
                       description="ID of parent SV of which this record is one part"),
                  dict(id='PARENT_SVTYPE', number=1, type='String',
                       description="type of parent SV of which this record is one part")
                  ]


def get_vcf_header_infos() -> list[dict]:
    vcf_header_infos = copy.deepcopy(VCF_HEADER_INFOS)
    for variant_set_class in VARIANT_SET_CLASSES:
        vcf_header_infos.extend(variant_set_class.get_vcf_header_infos())
    return vcf_header_infos

#############################################
