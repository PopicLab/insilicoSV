from constants import *
import utils
import random


class Structural_Variant():
    def __init__(self, sv_type, mode, length_ranges=None, source=None, target=None):
        '''
        Initializes SV's transformation and sets up its events, along with several other basic attributes like zygosity

        sv_type: Enum either specifying one of the prewritten classes or a Custom transformation, in which case source and target are required
        mode: flag to indicate whether we're in fixed or randomized mode
            ---> is there a better way to do this (subclass this and create an alternate constructor for the fixed case
                where a lot of this code isn't necessary)
        Params given in randomized mode:
        length_ranges: list containing tuple(s) (min_length, max_length) OR singleton int given if the SV is of known
            position/is being read in from a vcf of known/fixed variants
        source: tuple representing source sequence, optional
        target: tuple representing target sequence, optional
        '''
        self.type = sv_type
        if self.type != Variant_Type.Custom:
            self.source, self.target = SV_KEY[self.type]
            self.name = self.type.name
        else:
            self.source, self.target = tuple(source), tuple(target)
            self.name = str("".join(self.source)) + ">" + str("".join(self.target))

        # events like dispersions will appear as the same symbol, so it's important to add unique tags to differentiate them
        # user-generated source/target will also not group together duplication markings and their corresponding symbols together as it is inputted as string
        utils.validate_symbols(self.source, self.target)
        self.source_unique_char, self.target_unique_char = Structural_Variant.reformat_seq(
            self.source), Structural_Variant.reformat_seq(self.target)

        # initialize event classes
        self.start = None  # defines the space in which SV operates
        self.end = None
        self.req_space = None  # required space for SV, sum of event lengths
        self.source_events = []  # list of Event classes for every symbol in source sequence
        self.events_dict = dict()  # maps every unique symbol in source and target to an Event class
        # initialize_events sets the values of events_dict, source_dict, and req_space
        if mode == 'randomized':
            self.initialize_events(length_ranges)
        self.source_symbol_blocks = []
        self.target_symbol_blocks = []

        # specifies if sv is unable to be simulated due to random placement issues
        # will be turned on later
        self.active = False

        # 1 = homozygous, 0 = heterozygous
        self.ishomozygous = Zygosity.UNDEFINED

        # stores list of booleans specifying if SV will be applied to certain haplotype for assigned chromosome
        self.hap = [False, False]

    def __repr__(self):
        return "<SV transformation \"{}\" -> \"{}\" taking up {} non-dispersion spaces>".format(
            ''.join(self.source), ''.join(self.target),
            sum([event.length for event in self.source_events if not event.symbol.startswith(Symbols.DIS.value)]))

    @staticmethod
    def reformat_seq(transformation):
        # if dispersion events exist in transformation, tag on unique ids to make them distinct as they all are "_"
        # unique ids necessary to map symbol to event with one-to-one correspondence
        # groups together a duplication marking and its corresponding symbol
        # transformation: tuple, user inputted source and target
        unique_transform = []
        unique_id = 1
        for component in transformation:
            if component != Symbols.DIS.value and component != Symbols.DUP_MARKING.value:
                unique_transform.append(component)
            elif component == Symbols.DUP_MARKING.value:  # duplication event case, need to group together symbol and duplication marking
                unique_transform[-1] += Symbols.DUP_MARKING.value
            else:  # dispersion event case, component = dispersion
                unique_transform.append(component + str(unique_id))
                unique_id += 1
        return tuple(unique_transform)

    def initialize_events(self, lengths):
        '''
        Initializes event classes and creates a mapping of symbol to event

        lengths: list of tuples specifying min and max length for events within SV
        -> returns list of events in source sequence
        '''
        # collect all unique symbols present in both source and target sequences - include target as there may be insertions
        # note that the symbols represent events
        all_symbols = []
        for ele in self.source_unique_char + self.target_unique_char:
            # only append original symbols or dispersion events
            if len(ele) > 0 and (len(ele) == 1 or ele.startswith(Symbols.DIS.value)) and ele.upper() not in all_symbols:
                all_symbols.append(ele.upper())
        all_symbols.sort()  # user inputs symbol lengths in lexicographical order

        # symbols_dict: (key = symbol, value = (chosen length, length range))
        # determine length of events/symbols
        symbols_dict = dict()
        if len(lengths) > 1:  # values given by user represents custom ranges for each event symbol of variant in lexicographical order
            assert (len(lengths) == len(all_symbols)), \
                "Number of lengths entered does not match the number of symbols (remember foreign insertions and dispersions) present!"
            for idx, symbol in enumerate(all_symbols):
                symbols_dict[symbol] = (random.randint(lengths[idx][0], lengths[idx][1]), lengths[idx])

        elif len(lengths) == 1:  # value given by user represents length (same range) of each event within variant in lexicographical order
            for symbol in all_symbols:
                symbols_dict[symbol] = (random.randint(lengths[0][0], lengths[0][1]), lengths[0])

        else:
            raise Exception("Lengths parameter expects at least one tuple")
        # symbols_dict[Symbols.PLACEHOLDER] = (0, (0,0))

        # initialize event classes
        for idx, symbol in enumerate(all_symbols):
            # empty event - no source fragment yet
            event = Event(self, symbols_dict[symbol][0], symbols_dict[symbol][1], symbol)
            self.events_dict[symbol] = event

        for symbol in self.source_unique_char:
            self.source_events.append(self.events_dict[symbol])

        self.req_space = sum([event.length for event in self.source_events])
        # ** I don't think this return statement is used
        # return self.source_events

    def generate_blocks(self):
        '''
        Groups together source and target symbols between dispersion events (_)
        Tracks which block index the original symbols are in

        -> returns list of lists
        '''

        def find_blocks(transformation):
            # transformation: tuple of strings
            # -> returns list of lists of strings
            # Ex. ("A","B","_","C","D") -> [["A","B"], ["C","D"]]
            # Ex. ("A","B","_","_") -> [["A","B"],[],[]]
            blocks = [[]]
            for symbol in transformation:
                if not symbol.startswith(Symbols.DIS.value):
                    blocks[-1].append(symbol)
                else:
                    blocks.append([])
            return blocks

        def track_original_symbol(symbol_blocks):
            # finds which region/block the original symbol is in
            # if symbol is later found in another region, then translocation detected
            # blocks: list of lists of symbols
            for idx, block in enumerate(symbol_blocks):
                for symbol in block:
                    if len(symbol) == 1:  # means it's an original symbol
                        self.events_dict[symbol].original_block_idx = idx

        self.source_symbol_blocks = find_blocks(self.source_unique_char)
        self.target_symbol_blocks = find_blocks(self.target_unique_char)
        track_original_symbol(self.source_symbol_blocks)

        return self.target_symbol_blocks

    def change_fragment(self):
        '''
        Takes the mapping of symbols to events and the target sequence to construct a replacement sequence for the reference fragment
        '''
        decode_funcs = {"invert": lambda string: utils.complement(string[::-1]),
                        "identity": lambda string: string,
                        "complement": utils.complement}
        encoding = self.events_dict  # maps symbol like A or B to base pairs on reference

        # find all blocks of symbols between dispersion events
        # we will apply edits based on a block's start and end pos
        self.generate_blocks()  # blocks are the groupings of non-dispersion events

        changed_fragments = []
        assert (self.start is not None and self.end is not None), "Undefined SV start for {}".format(
            self)  # start & end should have been defined alongside event positions
        block_start = self.start  # describes SV's start position - applies for the first "block"
        curr_chr = self.start_chr

        for idx, block in enumerate(self.target_symbol_blocks):
            new_frag = ""
            for x, ele in enumerate(block):
                # used to find corresponding event from encoding, all keys in encoding are in uppercase
                upper_str = ele[0].upper()
                event = encoding[upper_str[0]]

                if any(c.islower() for c in ele):  # checks if lowercase symbols exist in ele, represents an inversion
                    new_frag += decode_funcs["invert"](event.source_frag)

                elif upper_str[0] in encoding:  # take original fragment, no changes
                    new_frag += event.source_frag

                elif ele.startswith(Symbols.DIS.value):  # DIS = dispersion event ("_")
                    raise Exception("Dispersion event detected within block: {}".format(self.target_symbol_blocks))
                else:
                    raise Exception("Symbol {} failed to fall in any cases".format(ele))

            # find dispersion event right after block to find position of next block
            assert curr_chr != None, "Unvalid chr detected for SV {} and events_dict {}".format(self, self.events_dict)
            if idx < len(self.target_symbol_blocks) - 1:
                dis_event = self.events_dict[Symbols.DIS.value + str(idx + 1)]  # find the nth dispersion event
                changed_fragments.append(
                    [curr_chr, block_start, dis_event.start, new_frag])  # record edits going by block
                block_start = dis_event.end  # move on to next block
                curr_chr = dis_event.source_chr
            else:
                changed_fragments.append([curr_chr, block_start, self.end, new_frag])

        self.changed_fragments = changed_fragments
        self.clean_event_storage()  # clean up unused storage - we do not need to store most source_frags anymore

        return changed_fragments

    def clean_event_storage(self):
        # remove source fragments from events to save space as they are no longer needed
        for event in self.events_dict.values():
            if event.symbol in self.source_unique_char:  # do not clean out insertion fragments as they'll need to be exported later
                event.source_frag = "Removed"


class Event():
    '''represents the symbols, also known as the "events," within a SV transformation'''

    def __init__(self, sv_parent, length, length_range, symbol):
        '''
        sv_parent: Structural Variant, event is always part of larger SV
        '''
        self.sv_parent = sv_parent
        self.length = length
        self.length_range = length_range
        self.symbol = symbol  # refers to symbol in SV's transformation
        self.source_chr = None
        self.source_frag = None
        self.start = None
        self.end = None

    def __repr__(self):
        return "<Event {}>".format({"length": self.length, "symbol": self.symbol, "start": self.start, "end": self.end,
                                    "source_chr": self.source_chr, "source_frag": self.source_frag})
