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
from insilicosv.utils import RegionFilter, OverlapMode, Locus, error_context, chk, TandemRepeatRegionFilter, if_not_none
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
        chk(isinstance(vset_config, dict), f'Each variant set must be a dict, {vset_config} provided',
            error_type='value')
        self.vset_config = copy.deepcopy(vset_config)
        self.config = config
        self.overlap_kinds = utils.as_list(self.vset_config.get('overlap_region_type', 'all'))
        self.overlap_ranges = []
        self.header = []

        self.vset_config['config_descr'] = ', '.join("%s: %s" % item for item in self.vset_config.items())
        self.novel_insertion_seqs = None

        # For predefined types given as custom in the config file, to keep track they have to be treated as custom
        self.input_type = None
        self.overlap_mode = None

        # For SNPs and INDELs overlap
        self.overlap_sv = False

        if 'overlap_mode' in self.vset_config:
            chk(isinstance(self.vset_config['overlap_mode'], str) or
                (isinstance(self.vset_config['overlap_mode'], list) and all(
                    isinstance(ovlp, str) for ovlp in self.vset_config['overlap_mode'])),
                f'Invalid overlap_mode, only a string or list of strings are accepted {vset_config}', error_type='syntax')
            try:
                self.overlap_mode = OverlapMode(self.vset_config['overlap_mode'])
            except ValueError:
                chk(False, f'Invalid overlap_mode in {vset_config}', error_type='value')

        self.vset_config['overlap_region_type'] = (tuple(utils.as_list(self.vset_config['overlap_region_type']))
                                                   if 'overlap_region_type' in self.vset_config else ('all',))
        if self.overlap_mode is not None and 'overlap_region_type' not in self.vset_config:
            logger.warning(f'overlap_region_type is not specified for {self} all provided ROIs will be considered')
        self.overlap_ranges = tuple(
            self.vset_config['overlap_region_length_range'] if 'overlap_region_length_range' in self.vset_config else
            [None, None])

        self.svtype = None
        if 'type' in self.vset_config:
            self.vset_config['type'] = self.vset_config['type'].replace(" ", "")
            if '->' in self.vset_config['type']:
                self.svtype = VariantType.CUSTOM
                grammar = self.vset_config['type'].split('->')
            else:
                self.svtype = VariantType(self.vset_config['type'])
                grammar = SV_KEY[self.svtype]

            self.source = grammar[0]
            chk(not Syntax.DIVERGENCE in self.source and not Syntax.MULTIPLE_COPIES in self.source,
                f'Operators ({Syntax.MULTIPLE_COPIES} and {Syntax.DIVERGENCE}) '
                f'are not allowed in the source, {vset_config}', error_type='syntax')

            self.target = grammar[1]
            chk(not Syntax.ANCHOR_START in self.target and not Syntax.ANCHOR_END in self.target,
                f'Anchors ({Syntax.ANCHOR_START} and {Syntax.ANCHOR_END}) '
                f'are not allowed in the target, If for instance, you want to constrain the insertion target of a dDUP, '
                f'write A_(). Error in {vset_config}', error_type='syntax')

        self.interchromosomal = self.vset_config.get('interchromosomal', False)
        chk(isinstance(self.interchromosomal, bool), 'interchromosomal if provided must be a boolean. '
                                                     'But %s was provided' % vset_config, error_type='syntax')
        self.interchromosomal_period = self.vset_config.get('interchromosomal_period', None)
        if self.interchromosomal and self.interchromosomal_period is None:
            self.interchromosomal_period = 0

        if self.interchromosomal_period is not None:
            self.interchromosomal = True

    @staticmethod
    def get_sampled_int_value(value, locals_dict=None):
        if isinstance(value, (int, float)) or value is None:
            return value
        elif isinstance(value, list):
            chk(len(value) == 2, f'Expected [min, max] pair: {value}', error_type='value')
            chk(isinstance(value[0], int) and isinstance(value[1], int) and
                int(value[0]) <= int(value[1]), f'Invalid [min, max] pair: {value}', error_type='value')
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
                chk(False, f'Error valuating the expression {value}: {exc}', error_type='value')

    @staticmethod
    def get_sampled_float_value(value):
        if isinstance(value, (int, float)):
            return value
        elif isinstance(value, list):
            chk(len(value) == 2, f'Expected [min, max] pair: {value}', error_type='value')
            chk(isinstance(value[0], float) and isinstance(value[1], float) and
                int(value[0]) <= int(value[1]), f'Invalid [min, max] pair: {value}', error_type='value')
            return random.uniform(value[0], value[1])
        chk(False, f'Error valuating the expression {value}', error_type='value')

    def grammar_to_variant_set(self, lhs_strs, rhs_strs, symbol_lengths, symbol_min_lengths, num_letters,
                               novel_insertion_seqs, n_copies_list, divergence_prob_list, replacement_seq=None, orig_seq=None,
                               vset_config=None,):
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
                breakend_interval_lengths.append(symbol_lengths[num_letters + len(lhs_dispersion)])
                breakend_interval_min_lengths.append(symbol_min_lengths[num_letters + len(lhs_dispersion)])
                lhs_dispersion.append(len(lhs))
            else:
                chk((len(lhs_str) == 1 and lhs_str[0].isupper()), f'Invalid LHS symbol {lhs_str} in the left hand side of {vset_config}')
                chk(len(symbol_lengths) >= len(lhs) - len(lhs_dispersion), f'Not enough length ranges for the left hand side in {vset_config}')
                breakend_interval_lengths.append(symbol_lengths[len(lhs) - len(lhs_dispersion)])
                breakend_interval_min_lengths.append(symbol_min_lengths[len(lhs) - len(lhs_dispersion)])
            lhs.append(Symbol(lhs_str))
        chk(((overlap_anchor_bounds[0] is None) and (overlap_anchor_bounds[1] is None)) or
            ((overlap_anchor_bounds[0] is not None) and (overlap_anchor_bounds[1] is not None)),
            f'Misformed anchor: {lhs_strs} in {vset_config}', error_type='syntax')
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
        # If a custom event specifies letters that are not modified and only use for placing constraints, we want ot keep track of them
        identities_to_add = {letter: True for letter in lhs if letter.name != Syntax.DISPERSION}
        # Parse the rhs to define the operations.
        for rhs_str in rhs_strs:
            if rhs_str == Syntax.DISPERSION:
                n_dispersions_rhs += 1
                chk(n_dispersions_rhs - 1 < len(lhs_dispersion),
                    f'The number of dispersions needs to be the same in the left and right hand sides. {lhs_strs} -> {rhs_strs}'
                    f' in {vset_config}',
                    error_type='syntax')
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
                                                            f'for each `{Syntax.MULTIPLE_COPIES}` symbol used. '
                                                            f'Error in {vset_config}', error_type='syntax')
                n_copies = self.get_sampled_int_value(n_copies_list[n_multiple_copies])
                chk(n_copies >= 1, f'The number of copies must be strictly positive, got {n_copies} for '
                                   f'the {n_multiple_copies + 1} \'{Syntax.MULTIPLE_COPIES}\' symbol {symbol}. Error in {vset_config}', error_type='syntax')
                n_multiple_copies += 1

            # Determine the divergence probability for each divergence symbol.
            divergence_prob = 0
            if Syntax.DIVERGENCE in rhs_str:
                divergence_prob = 1.

                if not replacement_seq:
                    # replacement_seq is used to provide the SNP when loading from vcf
                    chk(n_divergence_prob < len(divergence_prob_list), f'A divergence must be provided '
                                                                       f'for each \'{Syntax.DIVERGENCE}\' symbol used. '
                                                                       f'Error in {vset_config}', error_type='syntax')
                    divergence_prob = self.get_sampled_float_value(divergence_prob_list[n_divergence_prob])
                    chk(0 < divergence_prob <= 1,
                        f'The divergence probability must be between 0 (excluded) and 1 (included), got {divergence_prob} for '
                        f'the {n_divergence_prob + 1}-th \'{Syntax.DIVERGENCE}\' symbol {symbol}. Error in {vset_config}', error_type='syntax')

                if self.svtype and self.svtype != VariantType.SNP:
                    # A divergence can only be applied to a copied sequence
                    chk("".join(rhs_strs).upper().count(symbol.name) > 1,
                        f'{Syntax.DIVERGENCE} can only be applied to a duplicated sequence. '
                        f'Error in the {n_divergence_prob + 1}-th {Syntax.DIVERGENCE} of {vset_config}', error_type='syntax')
                    chk(divergence_prob < 1, f'The divergence probability of a divergence must be strictly less than 1. '
                                             f'Error in the {n_divergence_prob + 1}-th of {vset_config}', error_type='syntax')
                n_divergence_prob += 1

            transform = Transform(
                transform_type=transform_type,
                is_in_place=is_in_place,
                divergence_prob=divergence_prob,
                n_copies=n_copies,
                replacement_seq=replacement_seq,
                orig_seq=orig_seq
            )
            # We do not add operations for inplace identity transformations without divergence (do not affect the sequence).
            if (
                    transform_type != TransformType.IDENTITY or not is_in_place or divergence_prob > 0 or replacement_seq is not None
                    or n_copies > 1):
                # The letter has been involved in an operation and is not a placeholder.
                identities_to_add[symbol] = False
                operation = Operation(transform=transform, op_info={'SYMBOL': symbol.name})
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
                            # and in order of appearance in the rhs
                            chk(len(symbol_lengths) >= len(lhs) - len(lhs_dispersion) + len(novel_insertions),
                                f'Missing length_ranges for {symbol} in {lhs_strs}', error_type='syntax')
                            novel_insertion = None
                            if novel_insertion_seqs is not None:
                                chk(len(novel_insertion_seqs) > len(novel_insertions),
                                    f'If a novel_insertions field is provided, it must '
                                    f'contain the sequences for each novel insertion or None.', error_type='syntax')
                                if novel_insertion_seqs[len(novel_insertions)] is not None:
                                    novel_insertion = novel_insertion_seqs[len(novel_insertions)]
                            if novel_insertion is None:
                                novel_insertion = utils.generate_seq(
                                    length=symbol_lengths[len(lhs) - len(lhs_dispersion) + len(novel_insertions)])
                            novel_insertions[symbol] = novel_insertion
                        operation.novel_insertion_seq = novel_insertions[symbol]
                operations.append(operation)
        # Add the deletions
        for letter, delete in delete_letters.items():
            if delete:
                # The symbol is involved in a DEL, it is not a placeholder.
                identities_to_add[letter] = False
                source_index = lhs.index(letter)
                # The symbol was deleted, we add a DEL operation.
                transform = Transform(transform_type=TransformType.DEL, is_in_place=True)
                operation = Operation(transform=transform, op_info={'SYMBOL': letter.name})
                operation.source_breakend_region = BreakendRegion(source_index, source_index + 1)
                operations.append(operation)
        # Add the identities for place holder letters
        for letter, identity in identities_to_add.items():
            if identity:
                source_index = lhs.index(letter)
                # The symbol was deleted, we add a DEL operation.
                transform = Transform(transform_type=TransformType.IDENTITY, is_in_place=True)
                operation = Operation(transform=transform, op_info={'SYMBOL': letter.name})
                operation.source_breakend_region = BreakendRegion(source_index, source_index + 1)
                operations.append(operation)
        chk(len(lhs_dispersion) == n_dispersions_rhs, f'Dispersions cannot be altered, the same dispersions must '
                                                      f'be present on the left and right sides. lhs: {len(lhs_dispersion)} dispersions,'
                                                      f'rhs: {n_dispersions_rhs} dispersions.', error_type='syntax')
        # If overlap_anchor is not specified, infer "full SV" as the anchor if no dispersion
        if (overlap_anchor is None) and (self.overlap_mode is not None) and (not lhs_dispersion):
            overlap_anchor = BreakendRegion(0, len(lhs) - 1)
        return operations, overlap_anchor, lhs_dispersion, breakend_interval_lengths, breakend_interval_min_lengths

    # end: def grammar_to_variant_set


class SimulatedVariantSet(VariantSet):

    def __init__(self, vset_config, config):
        super().__init__(vset_config, config)

        chk(isinstance(vset_config.get('type'), str), f'Missing or bad variant type in {vset_config}', error_type='syntax')
        chk((isinstance(vset_config.get('number'), int) and vset_config['number'] >= 0),
            f'Missing or bad "number" of variants to generate in {vset_config}', error_type='value')

        self.preprocess_config()

    def preprocess_config(self):
        chk('divergence_prob' not in self.vset_config or isinstance(self.vset_config['divergence_prob'], (list, int, float)),
            f'divergence_prob must be a float or an int or a list of floats in ]0, 1] or a list of ranges. But, a '
                 f'%s was provided in %s' % (type(self.vset_config.get('divergence_prob', [])), self.vset_config), error_type='value')

        if ((self.svtype != VariantType.CUSTOM) and (Syntax.MULTIPLE_COPIES in ''.join(self.target))
                and ('n_copies' not in self.vset_config)):
            # Default the number of copies to 1 for predefined types with duplications
            self.vset_config['n_copies'] = [1]

        chk(isinstance('n_copies' not in self.vset_config or self.vset_config['n_copies'], (list, int)),
            f'The number of copies must be an integer or a list of integers or a list of ranges in {self.vset_config}',
            error_type='value')
        if 'n_copies' in self.vset_config and isinstance(self.vset_config['n_copies'], int):
            self.vset_config['n_copies'] = [self.vset_config['n_copies']]

        if self.svtype == VariantType.mCNV:
            chk(('n_copies' in self.vset_config) and (self.vset_config['n_copies'][0] > 1),
                f'n_copies has to be provided and be above 1 for a mCNV in {self.vset_config}', error_type='value')

    # end: def preprocess_config(self)

    @property
    def sv_type(self):
        return self.svtype

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
        if (self.config.get('homozygous_only', False) or (random.randint(0, 1) and not
           self.config.get('heterozygous_only', False))):
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
        return ((vset_config.get("type") is not None) and
                (vset_config.get("type") in [k.value for k in SV_KEY.keys()] or '->' in vset_config.get("type")))

    @override
    def preprocess_config(self):
        super().preprocess_config()

        for vset_config_key in self.vset_config:
            chk(vset_config_key in (
                'type', 'number',
                'length_ranges',
                'overlap_region_type', 'overlap_region_length_range',
                'overlap_mode',
                'blacklist_region_type',
                'divergence_prob',
                'n_copies',
                'novel_insertions',
                'interchromosomal',
                'interchromosomal_period',
                'config_descr',
                'allow_sv_overlap',
                'VSET'
            ), f'invalid SV config key {vset_config_key}', error_type='syntax')

        vset_cfg = self.vset_config
        rhs_strs_list: list[str] = []
        for c in self.target:
            if c in (Syntax.DIVERGENCE, Syntax.MULTIPLE_COPIES):
                chk(rhs_strs_list and rhs_strs_list[-1] and rhs_strs_list[-1][0].isalpha(),
                    f'{c} must modify a letter {rhs_strs_list}.', error_type='syntax')
                rhs_strs_list[-1] += c
            else:
                rhs_strs_list.append(c)

        self.target = rhs_strs_list

        # Anchor the whole SV if overlap mode is specified, no anchor is provided and there is no dispersion.
        if ((self.overlap_mode is not None) and (Syntax.DISPERSION not in self.source) and
                (Syntax.ANCHOR_START not in self.source)):
            self.source = tuple([Syntax.ANCHOR_START, *self.source, Syntax.ANCHOR_END])

        # Check if a type provided as grammar is a predefined type
        if self.svtype == VariantType.CUSTOM:
            lhs = tuple([letter for letter in self.source if letter not in [Syntax.ANCHOR_END, Syntax.ANCHOR_START]])
            rhs = tuple([letter if Syntax.DIVERGENCE not in letter and Syntax.MULTIPLE_COPIES not in letter
                         else letter[0] for letter in self.target])

            for key, grammar in SV_KEY.items():
                grammar = (grammar[0], tuple([letter if Syntax.MULTIPLE_COPIES not in letter and Syntax.DIVERGENCE not in letter
                                              else letter[0] for letter in grammar[1]]))
                # Test if the grammar or its symmetric match the grammar of the record
                if (grammar[0] == lhs and grammar[1] == rhs) or (
                        grammar[0] == lhs[::-1] and grammar[1] == rhs[::-1]):
                    # Distinguish between SNP/mCNV/DIVERGENCE/Identity
                    if (key == VariantType.mCNV) and (Syntax.MULTIPLE_COPIES not in self.target[0]): continue

                    if (key == VariantType.SNP and ((vset_cfg.get('length_ranges') not in (None, [[1, 1]])) or (Syntax.DIVERGENCE not in self.target[0])
                                                    or ('divergence_prob' in vset_cfg and vset_cfg['divergence_prob'] not in [[1], 1, 1., [1.]]))):
                            continue
                    # we found a match and update the types
                    self.svtype = key
                    self.input_type = VariantType.CUSTOM
                    break

        chk(self.overlap_mode != OverlapMode.CHROM or self.svtype in [VariantType.DEL, VariantType.DUP],
            'Only DEL and DUP SVs are allowed '
            'to have overlap_mode: chrom. Error in %s' % vset_cfg['config_descr'], error_type='syntax')

        chk('divergence_prob' not in vset_cfg or Syntax.DIVERGENCE not in self.target,
            f'\'{Syntax.DIVERGENCE}\' is not used but divergence_prob has been provided in {vset_cfg}', error_type='syntax')

        if self.svtype == VariantType.SNP:
            chk(vset_cfg.get('length_ranges') in (None, [[1, 1]]),
                f'length_ranges for SNP can only be [[1, 1]]. Error in %s' % vset_cfg['config_descr'], error_type='value')
            chk('divergence_prob' not in vset_cfg or vset_cfg['divergence_prob'] in [[1], 1],
                f'divergence prob for SNP can only be 1. Error in %s' % vset_cfg['config_descr'], error_type='value')
            vset_cfg['length_ranges'] = [[1, 1]]
            vset_cfg['divergence_prob'] = [1.0]
        elif vset_cfg['type'] == 'INDEL':
            # INDELS have to be of size <= 50
            chk(not vset_cfg.get('length_ranges', False) or (1 <= if_not_none(vset_cfg['length_ranges'][0][0], 1) <=
                                                          if_not_none(vset_cfg['length_ranges'][0][1], 50) <= 50 ),
                f'length_ranges for INDEL must be included in [0, 50]. Error in %s' % vset_cfg['config_descr'],
                error_type='value')
            chk(not vset_cfg.get('overlap_region_length_range', False) or (if_not_none(vset_cfg['overlap_region_length_range'][0], 0) <=
                                                          if_not_none(vset_cfg['overlap_region_length_range'][1], 50) <= 50 ),
                f'overlap_region_length_range for INDEL must be included in [0, 50]. Error in %s' % vset_cfg['config_descr'],
                error_type='value')

            # In case the max length_range was null, set it to 50 as it is the maximum for INDEL
            if vset_cfg.get('length_ranges', False) and not vset_cfg['length_ranges'][0][1]:
                vset_cfg['length_ranges'][0][1] = 50
                if not vset_cfg['length_ranges'][0][0]:
                    vset_cfg['length_ranges'][0][0] = 1

            if not vset_cfg.get('length_ranges'):
                vset_cfg['length_ranges'] = [[1, 50]]

            if vset_cfg.get('overlap_mode', False):
                if not vset_cfg.get('overlap_region_length_range', False):
                    vset_cfg['overlap_region_length_range'] = [[1, 50]]


        if self.svtype in [VariantType.SNP, VariantType.INDEL] or (self.svtype in [VariantType.INS, VariantType.DEL] and
                                                                   'length_ranges' in vset_cfg and
                                                                   vset_cfg['length_ranges'][0][1] and
                                                                   vset_cfg['length_ranges'][0][1] < 50):
            self.overlap_sv = vset_cfg.get('allow_sv_overlap', False)
        else:
            chk(not ('allow_sv_overlap' in vset_cfg),
                f'overlap_sv are only available for SNPs and INDELs, but %s was provided' %
                vset_cfg['config_descr'], error_type='type')

            chk('length_ranges' in vset_cfg or 'novel_insertions' in self.vset_config,
                f'Please specify length ranges in %s' % (vset_cfg['config_descr']), error_type='syntax')

            if 'length_ranges' in vset_cfg:
                chk(isinstance(vset_cfg['length_ranges'], list), f'length_ranges must be a list for %s' % vset_cfg['config_descr'],
                    error_type='syntax')
                for length_range in vset_cfg['length_ranges']:
                    chk(isinstance(length_range, str) or
                        (isinstance(length_range, list) and len(length_range) == 2 and
                         isinstance(length_range[0], (type(None), int, str)) and
                         isinstance(length_range[1], (type(None), int, str))),
                        f'invalid length_ranges. it must be a list of 2-tuples of str or int. '
                        f'Error in %s' % vset_cfg['config_descr'], error_type='value')

        if 'novel_insertions' in vset_cfg:
            try:
                with open(vset_cfg['novel_insertions'], 'r') as sequences:
                    self.novel_insertion_seqs = [line.rstrip() for line in sequences]
                    chk(all(bool(re.match('^[TCGA]+$', line)) for line in self.novel_insertion_seqs),
                        f'The file novel_insertions %s' % self.vset_config['novel_insertions'] +
                        f' contains invalid characters. It must be a list of sequences.', error_type='value')
            except:
                chk(False, f'novel_insertion file %s must be a readable ' % vset_cfg['novel_insertions'] +
                           f'file containing a sequence per line.', error_type='file not found')
        chk(isinstance(vset_cfg.get('type'), (type(None), list, str, tuple)),
            '%s must be a string or list of strings'.format(vset_cfg.get('type')), error_type='syntax')

        chk('interchromosomal_period' not in vset_cfg or isinstance(vset_cfg['interchromosomal_period'], (int, list)),
            'interchromosomal_period must be an int or a list of ints. '
            'Provided %s' % vset_cfg['config_descr'],
            error_type='syntax')
        interchromosomal_period = vset_cfg.get('interchromosomal_period', None)

        if isinstance(interchromosomal_period, list):
            chk(len(interchromosomal_period) == 2 and isinstance(interchromosomal_period[0], int) and
                isinstance(interchromosomal_period[1], int),
                'interchromosomal_period when provided as a range must contain two integers. '
                'Provided %s' % vset_cfg['config_descr'], error_type='syntax')

    def symmetrize(self, lhs_strs, rhs_strs, letter_ranges):
        # Enforce the symmetry of the predefined SVs with duplications or dispersions.
        if (("DUP" in self.svtype.name or "TRA" in self.svtype.name or
             "iDEL" in self.svtype.name)
                and (not self.interchromosomal)
                and (not self.input_type)
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
                letter_ranges = letter_ranges[::-1]
        chk(all(len(rhs) < 3 for rhs in rhs_strs),
            f'The operators + and * cannot be used at the same time {rhs_strs}.',
            error_type='syntax')
        return lhs_strs, rhs_strs, letter_ranges

    # Randomly pick distances withing the ranges defined for each symbol
    @staticmethod
    def pick_symbol_lengths(length_ranges, dispersion_ranges, letter_indexes, vset_config= None):
        ranges = length_ranges + dispersion_ranges
        remaining_symbols = [i for i in range(len(ranges))]
        chk(all(isinstance(length_range, str) or (isinstance(length_range, list) and (len(length_range) == 2))
                for length_range in ranges), f'length_ranges must be a list of [min, max] pairs {ranges} in %s' % vset_config['config_descr'], error_type='syntax')
        # Keep track of the dependencies between the symbols length ranges.
        # A dependency is the index of the letter the range is depending on (possibly a different letter for the min and max bounds)
        # and the offset to those letter lengths.
        dependencies = {}
        symbol_lengths = {}
        symbol_min_lengths = {}
        while remaining_symbols:
            idx = remaining_symbols.pop(0)
            min_range, max_range = ranges[idx]
            if (isinstance(min_range, int) or min_range is None) and (isinstance(max_range, int) or max_range is None):
                # Both bounds are independent to other letter lengths.
                def assign_length(min_range, max_range, is_dispersion):
                    if max_range is None:
                        length = None
                        if min_range is not None:
                            chk(is_dispersion, f'Only dispersions can have min but not max length {ranges} in %s' % vset_config['config_descr'], error_type='syntax')
                        min_length = min_range
                    else:
                        chk(min_range is not None, f'max_length given but not min_length {ranges} in %s' % vset_config['config_descr'], error_type='syntax')
                        chk(min_range <= max_range, f'max bound less than min bound {ranges} in %s' % vset_config['config_descr'], error_type='syntax')
                        chk(min_range >= 0, f'min length cannot be negative {ranges} in %s' % vset_config['config_descr'], error_type='value')
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
                        letters_bound = [(idx, char) for idx, char in enumerate(format_bound) if
                                         isinstance(char, str) and char.isalpha()]
                        chk(all(letter in letter_indexes for _, letter in letters_bound),
                            f'The length of a symbol depends on {letters_bound} '
                            f'one of which is not define in neither the source nor target in %s' % vset_config['config_descr'], error_type='value')
                        computed_letter = 0
                        for idx_in_dependency, letter in letters_bound:
                            # Gets the index of the letter in the list of length ranges.
                            index = letter_indexes[letter]
                            if index in symbol_lengths:
                                chk(symbol_lengths[index] is not None,
                                    f'A symbol length cannot depend on an unbounded symbol in %s' % vset_config['config_descr'], error_type='syntax')
                                ope = ''
                                if (idx_in_dependency > 0) and (isinstance(format_bound[idx_in_dependency - 1], int)
                                                                or format_bound[
                                                                    idx_in_dependency - 1] in letter_indexes):
                                    # The previous character was an integer, so a multiplication symbol was omitted.
                                    ope = '*'
                                format_bound[idx_in_dependency] = ope + str(symbol_lengths[index])
                                computed_letter += 1
                            else:
                                if not idx in dependencies:
                                    dependencies[idx] = []
                                dependencies[idx].append(index)
                                if index in dependencies:
                                    chk(not idx in dependencies[index],
                                        f'There is a cyclic dependency in the length definitions {letter} in %s' % vset_config['config_descr'],
                                        error_type='syntax')
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
                                chk(False, f'The length of a symbol depends on an invalid operation {idx} in %s' % vset_config['config_descr'],
                                    error_type='syntax')
                            chk(ranges[idx][pos] >= 0,
                                f'The operation defining the {idx}th symbol length gives a negative bound in %s' % vset_config['config_descr'],
                                error_type='value')
                remaining_symbols.append(idx)
        lengths = sorted(symbol_lengths.items(), key=lambda x: x[0])
        min_lengths = sorted(symbol_min_lengths.items(), key=lambda x: x[0])
        return [l[1] for l in lengths], [m_l[1] for m_l in min_lengths]

    @override
    def simulate_sv(self) -> SV:
        lhs_strs = self.source
        rhs_strs = self.target
        svtype = self.svtype
        if self.vset_config['type'] == 'INDEL':
            anchor = Syntax.ANCHOR_START in lhs_strs
            svtype = VariantType('DEL' if random.randint(0, 1) else 'INS')
            lhs_strs, rhs_strs = SV_KEY[svtype]
            if anchor:
                lhs_strs = tuple(Syntax.ANCHOR_START) + lhs_strs + tuple(Syntax.ANCHOR_END)
        length_ranges = self.vset_config['length_ranges'] if 'length_ranges' in self.vset_config else []
        letter_ranges = length_ranges
        dispersion_ranges = []
        lhs_strs_no_anchor = [letter for letter in lhs_strs if letter not in [Syntax.ANCHOR_END, Syntax.ANCHOR_START]]

        # Find the length ranges corresponding to dispersions and those corresponding to letters
        dispersions = [idx for idx, letter in enumerate(lhs_strs_no_anchor) if letter == Syntax.DISPERSION]
        if dispersions:
            if svtype != VariantType.CUSTOM and not self.input_type:
                # The SV was provided from a predefined type, the dispersion lengths are last
                letter_ranges = length_ranges[:-1]
                dispersion_ranges = [length_ranges[-1]]
            else:
                # The SV was provided from the grammar, the dispersion lengths are at the last position in the source
                letter_ranges = [length_range for idx, length_range in enumerate(length_ranges) if
                                 idx not in dispersions]
                dispersion_ranges = [length_range for idx, length_range in enumerate(length_ranges) if
                                     idx in dispersions]

        lhs_strs, rhs_strs, letter_ranges = self.symmetrize(lhs_strs, rhs_strs, letter_ranges)
        letters = [letter for letter in lhs_strs if
                   letter not in [Syntax.ANCHOR_END, Syntax.ANCHOR_START, Syntax.DISPERSION]]
        chk(len(letters) == len(set(letters)), f'Duplicate LHS symbol {letters} in {lhs_strs} for {self.vset_config}', error_type='syntax')

        # Add novel insertion letters only appearing in the rhs.
        for letter in rhs_strs:
            if letter[0].upper() not in letters + [Syntax.ANCHOR_END, Syntax.ANCHOR_START, Syntax.DISPERSION]:
                chk(letter[0].isupper(), 'A novel insertion letter has to be uppercase. But, %s was provided' % self.vset_config)
                letters.append(letter[0].upper())
        chk(len(length_ranges) == len(letters) + len(dispersions),
            f'Mismatched length ranges, expected {len(letters) + len(dispersions)} provided {len(length_ranges)} for '
            f'{self.vset_config}', error_type='syntax')
        letter_indexes = {letter: index for index, letter in enumerate(letters)}

        # Compute the lengths of the different symbols and dispersions from the length ranges.
        symbol_lengths, symbol_min_lengths = self.pick_symbol_lengths(letter_ranges, dispersion_ranges, letter_indexes, self.vset_config)

        novel_insertion_seqs = self.novel_insertion_seqs
        n_copies_list = self.vset_config.get('n_copies', [])

        divergence_prob_list = self.vset_config.get('divergence_prob', [])
        if not isinstance(divergence_prob_list, list):
            divergence_prob_list = [divergence_prob_list]

        # Build the different operations and anchor, determine the breakends and the distance between them.
        (operations, anchor, dispersions, breakend_interval_lengths,
         breakend_interval_min_lengths) = self.grammar_to_variant_set(lhs_strs, rhs_strs, symbol_lengths,
                                                                      symbol_min_lengths, len(letters),
                                                                      novel_insertion_seqs, n_copies_list,
                                                                      divergence_prob_list, vset_config=self.vset_config)

        #
        # construct the SV object
        #
        info = self.construct_info(lhs_strs, rhs_strs)
        roi_filter = self.get_roi_filter()

        interchromosomal_period = self.get_sampled_int_value(self.interchromosomal_period)

        return BaseSV(sv_id=self.make_sv_id(),
                      breakend_interval_lengths=breakend_interval_lengths,
                      breakend_interval_min_lengths=breakend_interval_min_lengths,
                      interchromosomal_period=interchromosomal_period,
                      operations=operations,
                      anchor=anchor,
                      dispersions=dispersions,
                      overlap_mode=self.overlap_mode,
                      roi_filter=roi_filter,
                      blacklist_filter=self.get_blacklist_filter(),
                      fixed_placement=None,
                      info=info,
                      genotype=self.pick_genotype(),
                      allow_sv_overlap=self.vset_config.get('allow_sv_overlap', False),
                      config_descr=self.vset_config['config_descr'])

    def construct_info(self, lhs_strs, rhs_strs):
        sv_type_str = self.svtype.name
        source_str = ''.join(lhs_strs)
        target_str = ''.join(rhs_strs)
        grammar = f'{source_str}->{target_str}'
        return {'OP_TYPE': sv_type_str, 'GRAMMAR': grammar}


# end: class FromGrammarVariantSetMaker

class TandemRepeatVariantSet(SimulatedVariantSet):

    @override
    @classmethod
    def can_make_from(cls, vset_config):
        return vset_config.get("type") in ('trEXP', 'trCON')

    def __init__(self, vset_config, config):
        super().__init__(vset_config, config)
        chk('repeat_count_change_range' in self.vset_config,
            'Please specify repeat_count_change_range for tandem repeat SVs.', error_type='syntax')
        chk(utils.is_valid_closed_int_range(self.vset_config['repeat_count_change_range']),
            'Invalid repeat_count_change_range', error_type='value')

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
        info = dict(SVTYPE=self.svtype.value, TR_CHANGE=repeat_count_change)
        overlap_region_type = (tuple(utils.as_list(self.vset_config['overlap_region_type']))
                               if 'overlap_region_type' in self.vset_config else 'all')
        anchor = BreakendRegion(0, 1)
        self.overlap_mode = OverlapMode.EXACT

        chk(self.overlap_mode in [None, OverlapMode.EXACT], 'overlap_mode must be "exact" for Tandem Repeat',
            error_type='syntax')
        if self.svtype == VariantType.trEXP:
            breakend_interval_lengths = [None]
            operations = [Operation(transform=Transform(transform_type=TransformType.IDENTITY,
                                                        is_in_place=False,
                                                        n_copies=repeat_count_change),
                                    source_breakend_region=BreakendRegion(start_breakend=0, end_breakend=1),
                                    target_insertion_breakend=0,
                                    target_insertion_order=(0,),
                                    op_info={'SYMBOL': 'A'})]
            # We only need the repeat motif to be present once.
            roi_filter = TandemRepeatRegionFilter(min_num_repeats=1,
                                                  region_kinds=overlap_region_type)
            return TandemRepeatExpansionContractionSV(
                sv_id=self.make_sv_id(),
                breakend_interval_lengths=breakend_interval_lengths,
                breakend_interval_min_lengths=[None] * len(breakend_interval_lengths),
                interchromosomal_period=None,
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
        elif self.svtype == VariantType.trCON:
            breakend_interval_lengths = [None]
            operations = [Operation(transform=Transform(transform_type=TransformType.DEL,
                                                        is_in_place=True, n_copies=1),
                                    source_breakend_region=BreakendRegion(0, 1),
                                    op_info={'SYMBOL': 'A'})]

            # ensure there are enough existing repeats to delete.
            roi_filter = TandemRepeatRegionFilter(min_num_repeats=repeat_count_change,
                                                  region_kinds=overlap_region_type)

            return TandemRepeatExpansionContractionSV(
                sv_id=self.make_sv_id(),
                breakend_interval_lengths=breakend_interval_lengths,
                breakend_interval_min_lengths=[None],
                interchromosomal_period=False,
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
    op_types = ['IDENTITY', 'COPY-PASTE', 'CUT-PASTE', 'COPYinv-PASTE', 'CUTinv-PASTE', 'NA', 'CUT']
    can_import_types = [v.value for v in VariantType] + op_types

    next_import_id: ClassVar[int] = 0

    @override
    @classmethod
    def can_make_from(cls, vset_config):
        return 'import' in vset_config

    def __init__(self, vset_config, config):
        super().__init__(vset_config, config)
        chk(utils.is_readable_file(vset_config['import']),
            '{path} vcf must name a readable file'.format(path=vset_config['import']), error_type='file not found')
        chk(set(vset_config.keys()) <= {'import', 'VSET'}, f'invalid config key in {vset_config}',
            error_type='syntax')

        self.import_id = ImportedVariantSet.make_import_id()
        with FastaFile(config['reference']) as reference:
            self.chrom_lengths = {chrom: chrom_length
                                  for chrom, chrom_length in zip(reference.references,
                                                                 reference.lengths)}

    @staticmethod
    def make_import_id():
        import_id = str(ImportedVariantSet.next_import_id)
        ImportedVariantSet.next_import_id += 1
        return import_id

    @override
    def make_variant_set(self):
        recs = defaultdict(list)
        num_simple_sv = 0
        with closing(VariantFile(self.vset_config['import'])) as vcf:
            self.header = vcf.header
            for vcf_rec in vcf.fetch():
                with error_context(vcf_rec):
                    chk(vcf_rec.chrom in self.chrom_lengths, 'An imported SV belong to a chromosome not'
                                                             'represented in the reference file.', error_type='value')
                    vcf_info = dict(vcf_rec.info)
                    if 'SVID' in vcf_info:
                        # Use the parent ID
                        recs[vcf_info['SVID']].append(vcf_rec)
                    elif vcf_rec.id:
                        recs[vcf_rec.id].append(vcf_rec)
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
        chk(vcf_rec.chrom in self.chrom_lengths, f'Unknown contig {vcf_rec.chrom}  in {vcf_rec}', error_type='value')
        chk(3 >= len(vcf_rec.alleles) >= 2, f'Only diploids are supported. But, {vcf_rec} was provided.', error_type='value')
        chk(len(vcf_rec.samples) <= 1, f'Can only import VCFs with one sample for {vcf_rec}', error_type='value')
        parsed_info = {}
        if vcf_rec.samples:
            sample = vcf_rec.samples[0]
            parsed_info['GENOTYPE'] = (bool(sample['GT'][0]), bool(sample['GT'][1]))
            if not sum(parsed_info['GENOTYPE']):
                logger.warning(f'The genotype provided in {vcf_rec} is invalid. A genotype will be randomly defined.')
                del parsed_info['GENOTYPE']
        if not 'GENOTYPE' in parsed_info:
            parsed_info['GENOTYPE'] = random.choice([(True, True), (True, False), (False, True)])
        vcf_info = dict(vcf_rec.info)

        parsed_info['ALLOW_SV_OVERLAP'] = vcf_info.get('ALLOW_SV_OVERLAP', False)

        if set(''.join(vcf_rec.alleles).upper().replace(' ', '')) <= set('TCGA'):
            if len(vcf_rec.alleles[0]) == 1 and 1 <= len(vcf_rec.alleles[1]) <= 2:
                if not vcf_info.get('SVTYPE'):
                    # SNP
                    vcf_info['OP_TYPE'] = 'SNP'
                    vcf_info['SVTYPE'] = 'SNP'

        chk('OP_TYPE' in vcf_info or 'SVTYPE' in vcf_info,
            f'Need an SVTYPE or OP_TYPE to import from vcf records for {vcf_rec}', error_type='syntax')
        rec_type_str = vcf_info.get('OP_TYPE', 'NA')

        chk(rec_type_str in {variant_type for variant_type in self.can_import_types},
            f'Currently only the following VCF types are supported: {self.can_import_types} but {rec_type_str} was provided',
            error_type='syntax')

        if rec_type_str not in self.op_types:
            parsed_info['OP_TYPE'] = VariantType(rec_type_str)
        else:
            parsed_info['OP_TYPE'] = rec_type_str

        chk(vcf_rec.start is not None, f'The Start position must not be None for {vcf_rec}', error_type='value')
        rec_start = Locus(chrom=vcf_rec.chrom, pos=vcf_rec.start)
        parsed_info['START'] = rec_start
        rec_end = Locus(chrom=vcf_rec.chrom, pos=vcf_rec.stop)

        if 'SVLEN' in vcf_info:
            rec_len = vcf_info['SVLEN']
            if isinstance(rec_len, (tuple, list)):
                chk(len(rec_len) == 1, f'Wrong format for the field SVLEN, only integers are supported. The record provided '
                                       f'was {vcf_rec}', error_type='syntax')
                rec_len = int(rec_len[0])
        else:
            rec_len = rec_end.pos - rec_start.pos
        parsed_info['END'] = rec_end
        parsed_info['SVLEN'] = rec_len

        if 'GRAMMAR' in vcf_info:
            parsed_info['GRAMMAR'] = vcf_info['GRAMMAR']

        if not (rec_len <= 50 and (rec_type_str == 'INS' or rec_type_str == 'INV' or vcf_info['SVTYPE'] == 'SNP' or
                              ('SVTYPE' in vcf_info and (vcf_info['SVTYPE'] == 'INS' or vcf_info['SVTYPE'] == 'INV')))):
            chk(not parsed_info['ALLOW_SV_OVERLAP'], f'ALLOW_SV_OVERLAP only allowed for SNPs or INDELs. But, {vcf_rec} was provided.')

        rec_target = None
        is_interchromosomal = False
        if 'TARGET' in vcf_info:
            chk(isinstance(vcf_info.get('TARGET'), int), f'TARGET has to be an int representing a position {vcf_info}',
                error_type='value')
            rec_target = Locus(chrom=vcf_info.get('TARGET_CHROM', vcf_rec.chrom), pos=vcf_info['TARGET'])
            if rec_target.chrom != vcf_rec.chrom:
                is_interchromosomal = True
            elif rec_end is not None:
                # The target has to be outside of the source region
                chk(rec_target.chrom != rec_start.chrom or rec_target.pos >= rec_end.pos or rec_target.pos <= rec_start.pos,
                    f"Invalid dispersion target {vcf_rec}", error_type='value')
        parsed_info['TARGET'] = rec_target
        parsed_info['INTERCHROMOSOMAL'] = is_interchromosomal

        novel_insertion_seq = None
        if 'INSSEQ' in vcf_info:
            chk(isinstance(vcf_info['INSSEQ'], str) and all(bp in 'TCGA' for bp in vcf_info['INSSEQ']),
                'INSSEQ has to be a valid sequence, {} was provided'.format(vcf_info['INSSEQ']), error_type='value')
            novel_insertion_seq = [vcf_info['INSSEQ']]
        elif parsed_info['OP_TYPE'] == VariantType.INS:
            # It is a novel sequence insertion but the sequence to insert isn't provided, we generate a random sequence
            logger.warning(
                'SV containing a novel sequence insertion but the sequence was not provided, a random sequence is used')
            novel_insertion_seq = [utils.generate_seq(length=rec_len)]

        if 'DIVERGENCE_PROB' in vcf_info and parsed_info['OP_TYPE'] != VariantType.SNP:
            logger.warning('InsilicoSV does not support importing specific divergence sequences for now,'
                           'a new mutated sequence will be created.')

        parsed_info['INSSEQ'] = novel_insertion_seq

        parsed_info['NCOPIES'] = vcf_info.get('NCOPIES', 1)
        parsed_info['INSORD'] = vcf_info.get('INSORD', None)
        parsed_info['DIVERGENCE_PROB'] = vcf_info.get('DIVERGENCE_PROB', [0])
        parsed_info['ALT'] = None
        parsed_info['REF'] = None
        parsed_info['SVTYPE'] = vcf_info.get('SVTYPE', 'Custom')

        if parsed_info['OP_TYPE'] == VariantType.SNP:
            parsed_info['DIVERGENCE_PROB'] = [1.0]
            if vcf_rec.alts[0] != '<SNP>':
                parsed_info['ALT'] = vcf_rec.alts
                chk(1 <= len(vcf_rec.alts) <= 2, f'Error in the ALT field format {vcf_rec}')
                parsed_info['ALT'] = [vcf_rec.alts[hap_index] if parsed_info['GENOTYPE'][hap_index] else None for hap_index in [0, -1]]
            if vcf_rec.ref != 'N':
                parsed_info['REF'] = vcf_rec.ref[0]
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
        parent_info = {'SVID': parent_id}
        positions_per_rec = []
        current_insord = 0

        # Extract the info of all the records of the SV to determine the breakends and placement.
        for vcf_rec in sv_recs:
            print(vcf_rec)
            parsed_info, additional_info = self.parse_vcf_rec_info(vcf_rec)

            if genotype is None:
                genotype = parsed_info['GENOTYPE']
            print(genotype, parsed_info)
            assert genotype == parsed_info['GENOTYPE']

            if len(sv_recs) > 1:
                parent_info['OP_TYPE'] = parsed_info['SVTYPE']
            else:
                if parsed_info['SVTYPE'] not in [VariantType.SNP, VariantType.CUSTOM]:
                    parsed_info['OP_TYPE'] = VariantType(parsed_info['SVTYPE'])
                parent_info['OP_TYPE'] = parsed_info['SVTYPE']

            parent_info['GRAMMAR'] = parsed_info.get('GRAMMAR', '')
            if 'GRAMMAR' not in parsed_info and parsed_info['SVTYPE'] != 'Custom':
                lhs_strs, rhs_strs = SV_KEY[VariantType(parsed_info['SVTYPE'])]
                parsed_info['GRAMMAR'] = ''.join(lhs_strs) + '->' + ''.join(rhs_strs)

            target_left = False
            source_regions.append([parsed_info['START'], parsed_info['END']])
            placement = [parsed_info['START'], parsed_info['END']]
            # Gather the length information to parse the grammar and get the operations associated to the record
            symbol_lengths = [parsed_info['END'].pos - parsed_info['START'].pos]
            symbol_min_lengths = [None]
            insord = (0,)
            if parsed_info['TARGET'] is not None:
                # We check the relative position of the target compared to start and end
                if parsed_info['END'].chrom == parsed_info['TARGET'].chrom and parsed_info['TARGET'] <= parsed_info['START']:
                    placement = [parsed_info['TARGET'], parsed_info['START'], parsed_info['END']]
                    target_left = True
                else:
                    placement = [parsed_info['START'], parsed_info['END'], parsed_info['TARGET']]
                    chk((parsed_info['END'].chrom != parsed_info['TARGET']) or (parsed_info['TARGET'] >= parsed_info['END']),
                        f'The position of the target has to be outside of the source region,'
                        f'{vcf_rec} has a target between the start and end.', error_type='value')
                if not parsed_info['INTERCHROMOSOMAL']:
                    # The target is an intrachromosomal dispersion
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
                    chk(parsed_info['INSORD'] is not None,
                        'Two sequences inserted at the locus {}, the INSORD field must be provided for these records to prevent any ambiguity.'.format(
                            parsed_info['INSORD']), error_type='value')
                insord = (parsed_info['INSORD'],) if parsed_info['INSORD'] is not None else (current_insord,)
                current_insord = insord[0] + 1

            insseq = parsed_info.get('INSSEQ', None)
            if insseq:
                insseq = insseq[0]

            targets.append(parsed_info['TARGET'])
            positions_per_rec.append(placement)
            operations = []
            is_in_place = parsed_info['TARGET'] is None
            # The SVTYPE of a record must be a predefined SV type or the identity operation.
            if ('GRAMMAR' in parsed_info and len(sv_recs) == 1) or (isinstance(parsed_info['OP_TYPE'], VariantType) and (
                    parsed_info['OP_TYPE'] not in [VariantType.DEL, VariantType.INV, VariantType.CUSTOM])):
                if 'GRAMMAR' in parsed_info:
                    chk(len(parsed_info['GRAMMAR'].split('->')) == 2, f'Unsupported GRAMMAR format {vcf_rec}.')
                    lhs_strs, rhs_strs = parsed_info['GRAMMAR'].split('->')
                else:
                    lhs_strs, rhs_strs = SV_KEY[parsed_info['OP_TYPE']]
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
                operations, _, _, _, _ = self.grammar_to_variant_set(lhs_strs, rhs_strs_list, symbol_lengths,
                                                                     symbol_min_lengths, 1,
                                                                     insseq,
                                                                     [parsed_info['NCOPIES']],
                                                                     parsed_info['DIVERGENCE_PROB'],
                                                                     replacement_seq=parsed_info['ALT'],
                                                                     orig_seq=parsed_info['REF'],
                                                                     vset_config=vcf_rec)
                for operation in operations:
                    operation.op_info = additional_info
                    operation.target_insertion_order = insord
            else:
                # The record is an atomic operation of a complex SV (CUT, COPY, INV...)
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

                op_attributes = []
                op_type = parsed_info['OP_TYPE']
                chk(op_type != 'NA', f'A custom type has been provided without OP_TYPE {vcf_rec}')
                if isinstance(parsed_info['OP_TYPE'], VariantType):
                    op_type = parsed_info['OP_TYPE'].name
                if ('inv' in op_type) or ('INV' in op_type):
                    op_attributes.append(('INV', is_in_place, target_breakend))
                elif ('PASTE' in op_type) or ('IDENTITY' in op_type):
                    op_attributes.append(('IDENTITY', is_in_place, target_breakend))
                if ('DEL' in op_type) or ('CUT' in op_type):
                    op_attributes.append(('DEL', True, None))
                for op_type, op_is_in_place, op_target in op_attributes:
                    transform = Transform(TransformType[op_type], is_in_place=op_is_in_place,
                                          n_copies=parsed_info['NCOPIES'],
                                          divergence_prob=parsed_info['DIVERGENCE_PROB'][0],
                                          replacement_seq=parsed_info['ALT'])
                    operations.append(Operation(transform, source_breakend_region=source_region,
                                                novel_insertion_seq=insseq,
                                                target_insertion_breakend=op_target,
                                                target_insertion_order=insord, op_info=additional_info))
            sv_operations.append(operations)
        # If we have more than one record we unify the breakends and their positions through the different operations
        if len(sv_recs) > 1:
            # Combine the information from the different records to determine the sets of breakends of the whole SV
            # Remove duplicate positions and insert the target positions in order to find the positions of all the breakends
            placements = self.remove_duplicates(
                [source[region_idx] for source in source_regions for region_idx in [0, 1]] +
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
                    ope_target = breakends[
                        operation.target_insertion_breakend] if operation.target_insertion_breakend is not None else None
                    start_breakend = breakends[operation.source_breakend_region.start_breakend]
                    end_breakend = breakends[operation.source_breakend_region.end_breakend]
                    operation.source_breakend_region = BreakendRegion(start_breakend=start_breakend,
                                                                      end_breakend=end_breakend)
                    operation.target_insertion_breakend = ope_target
            sv_operations = [sv_operation for operation_list in sv_operations for sv_operation in operation_list]
        else:
            # There is a single record
            sv_operations = sv_operations[0]
            placements = positions_per_rec[0]

        sv_id = 'Imported_' + self.import_id + '_' + str(parent_id)
        placement_dict = {breakend: locus for breakend, locus in enumerate(placements)}

        return BaseSV(sv_id=sv_id,
                      breakend_interval_lengths=[None] * (len(placements) - 1),
                      # Positions are known, the lengths are not needed
                      breakend_interval_min_lengths=[None] * (len(placements) - 1),
                      interchromosomal_period=None,  # The target chromosome is already known from the fixed_placement
                      operations=sv_operations,
                      fixed_placement=placement_dict,
                      overlap_mode=None,
                      anchor=None,
                      roi_filter=None,
                      blacklist_filter=None,
                      info=parent_info,
                      genotype=genotype,
                      allow_sv_overlap=parsed_info['ALLOW_SV_OVERLAP'],
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
    for variant_set_class in VARIANT_SET_CLASSES:
        if variant_set_class.can_make_from(vset_config):
            variant_set = variant_set_class(vset_config, config)
            return variant_set.make_variant_set(), variant_set.overlap_ranges, variant_set.overlap_kinds, variant_set.overlap_mode, variant_set.header
    chk(False, f"The format of the config or the sv_type is not supported {vset_config}")


VCF_HEADER_INFOS = [
    dict(id='END', number=1, type='Integer', description="End position of the variant described in this record"),
    dict(id='INSSEQ', number=1, type='String', description="Novel insertion sequence"),
    dict(id='INSORD', number=1, type='Integer',
         description="Insertion order for insertions at given position"),
    dict(id='VSET', number=1, type='Integer',
         description="Variant set number (numbered from 0) from which this variant was created"),
    dict(id='OP_TYPE', number=1, type='String',
         description="Type of operation performed for the corresponding record"),
    dict(id='GRAMMAR', number=1, type='String',
         description="Grammar of the structural variant"),
    dict(id='SYMBOL', number=1, type='String',
         description="Symbol the record considers"),
    dict(id='SVLEN', number=1, type='Integer',
         description="Length of structural variant"),
    dict(id='NCOPIES', number=1, type='Integer',
         description="Number of sequence copies to insert at target"),
    dict(id='DIVERGENCE_PROB', number=1, type='Float',
         description="Mutation probability for each nucleotide of a duplicated sequence."),
    dict(id='TARGET', number=1, type='Integer',
         description="Target location for a dispersed duplication or translocation"),
    dict(id='TARGET_CHROM', number=1, type='String',
         description="Target chromosome for a dispersed duplication or translocation"),
    dict(id='OVLP', number=1, type='String',
         description="Type of ROI on which the SV component was placed"),
    dict(id='OVLP_TARGET', number=1, type='String',
         description="Type of ROI on which the insertion target of an SV component was placed"),
    dict(id='OVLP_TYPE', number=1, type='String',
         description="Type of overlap with the ROI"),
    dict(id='ALLOW_SV_OVERLAP', number=1, type='String',
         description="If this record was allowed to overlap with other SVs."),
    dict(id='SVID', number=1, type='String',
         description="ID of parent SV of which this record is one part"),
    dict(id='SVTYPE', number=1, type='String',
         description="type of parent SV of which this record is one part")
    ]


def get_vcf_header_infos() -> list[dict]:
    vcf_header_infos = copy.deepcopy(VCF_HEADER_INFOS)
    for variant_set_class in VARIANT_SET_CLASSES:
        vcf_header_infos.extend(variant_set_class.get_vcf_header_infos())
    return vcf_header_infos

#############################################