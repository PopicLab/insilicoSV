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
        self.changed_fragments = []  # list recording new fragments to be placed in the output ref
        # initialize_events sets the values of events_dict, source_dict, and req_space
        if mode == 'randomized':
            self.initialize_events(length_ranges)

        sv_blocks = Blocks(self)
        # debug
        print(f'==== sv_blocks object ====\n{sv_blocks}')
        self.source_symbol_blocks = sv_blocks.source_blocks
        self.target_symbol_blocks = sv_blocks.target_blocks

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
            if component != Symbols.DIS.value and component != Symbols.DUP_MARKING.value and component != Symbols.DIV.value:
                unique_transform.append(component)
            elif component == Symbols.DUP_MARKING.value:  # duplication event case, need to group together symbol and duplication marking
                unique_transform[-1] += Symbols.DUP_MARKING.value
            elif component == Symbols.DIV.value: # divergence event case, want to keep track of interval needing modification
                unique_transform[-1] += Symbols.DIV.value
            else:  # dispersion event case, component = dispersion
                unique_transform.append(component + str(unique_id))
                unique_id += 1
        return tuple(unique_transform)

    @staticmethod
    def get_event_frag(event, symbol):
        # helper fn to get the ref frag for a given subevent
        # event: source event from events_dict
        # symbol: target symbol
        decode_funcs = {"invert": lambda string: utils.complement(string[::-1]),
                        "identity": lambda string: string,
                        "complement": utils.complement,
                        "diverge": lambda string: utils.divergence(string)}
        if any(c.islower() for c in symbol):
            return decode_funcs["invert"](event.source_frag)
        elif symbol[-1] == '*':  # checks if the element ends in an *, representing a divergent duplicate
            return decode_funcs["diverge"](event.source_frag)
        else:  # take original fragment, no changes
            return event.source_frag

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

        # initialize event classes
        for idx, symbol in enumerate(all_symbols):
            # empty event - no source fragment yet
            event = Event(self, symbols_dict[symbol][0], symbols_dict[symbol][1], symbol)
            self.events_dict[symbol] = event

        for symbol in self.source_unique_char:
            self.source_events.append(self.events_dict[symbol])

        self.req_space = sum([event.length for event in self.source_events])

        # debug
        print('END OF INITIALIZE_EVENTS')
        print(f'self.events_dict = {self.events_dict}')

    def assign_locations(self, start_pos):
        """
        assign events start and end positions (the single-sv-level version of choose_rand_pos originally)
        --> to be performed once source and target Blocks objects are populated and in the right order
        --> (so the process will be to choose a start position for the SV and assign all the successive positions accordingly)
        ** NB: will just apply this to target_blocks (the source representation doesn't appear to be necessary beyond
        being queried for change_fragment())
        ** usage: going to keep choose_rand_pos() in the simulator object and have that perform the drawing/checking of
        a valid start position (since that needs access to the max_tries setting from the config), but going to invoke
        this method to actually modify the target_blocks' events objects to set the start/end positions
        """
        current_pos = start_pos
        for block in self.target_symbol_blocks:
            for ev in block:
                if ev.symbol.startswith(Symbols.DIS.value):
                    current_pos += ev.length
                else:
                    ev.start = current_pos
                    ev.end = current_pos + ev.length
                    source_event = self.events_dict[ev.symbol[0].upper()]
                    ev.source_frag = self.get_event_frag(source_event, ev.symbol)
                    ev.source_chr = self.start_chr
                    current_pos = ev.end
        # debug
        print(f'===LOCATIONS ASSIGNED===\ntarget_symbol_blocks: {self.target_symbol_blocks}')

    def change_fragment(self):
        '''
        Takes the mapping of symbols to events and the target sequence to construct a replacement sequence for the reference fragment
        '''
        changed_fragments = []
        assert (self.start is not None and self.end is not None), "Undefined SV start for {}".format(
            self)  # start & end should have been defined alongside event positions
        block_start = None
        block_end = None

        # debug
        print('===CHANGE_FRAGMENT===')
        for block in self.target_symbol_blocks:
            new_frag = ''
            if block[0].symbol.startswith(Symbols.DIS.value):
                continue
            for i in range(len(block)):
                ev = block[i]
                new_frag += ev.source_frag
                if i == 0:
                    block_start = ev.start
                if i == len(block) - 1:
                    block_end = ev.end
            changed_fragments.append([self.start_chr, block_start, block_end, new_frag])
            print(f'new change fragment : {changed_fragments[-1]}')

            # # find dispersion event right after block to find position of next block
            # assert curr_chr != None, "Unvalid chr detected for SV {} and events_dict {}".format(self, self.events_dict)
            # if idx < len(self.target_symbol_blocks) - 1:
            #     dis_event = self.events_dict[Symbols.DIS.value + str(idx + 1)]  # find the nth dispersion event
            #     # debug
            #     # ** Here dis_event.start is being used as "block_end" but in the case of a flipped dispersion, these
            #     # ** won't be equal -- where can we establish the actual block_end?
            #     print(f'appending to changed fragments (idx={idx}):\n{str([curr_chr, block_start, block_end, new_frag])}')
            #     changed_fragments.append(
            #         # [curr_chr, block_start, dis_event.start, new_frag])  # record edits going by block
            #         [curr_chr, block_start, block_end, new_frag])  # <- block_end being the end of the last symbol in the block
            #     block_start = dis_event.end  # move on to next block
            #     print(f'block_start changed to {block_start}')
            #     curr_chr = dis_event.source_chr
            # else:
            #     # debug
            #     print(f'appending to changed fragments (idx={idx}):\n{str([curr_chr, block_start, block_start, new_frag])}')
            #     # ** for flipped dispersion self.end is in the wrong place (it should match with block_start)
            #     # ** --> it gets set to sv.start + sv.req_space in choose_rand_pos()
            #     changed_fragments.append([curr_chr, block_start, self.end, new_frag])
            #     # ** this is what we want for dDUPs/INV_dDUPs but not for events without dispersions
            #     # changed_fragments.append([curr_chr, block_start, block_start, new_frag])

        self.changed_fragments = changed_fragments
        self.clean_event_storage()  # clean up unused storage - we do not need to store most source_frags anymore

    def clean_event_storage(self):
        # remove source fragments from events to save space as they are no longer needed
        for event in self.events_dict.values():
            if event.symbol in self.source_unique_char:  # do not clean out insertion fragments as they'll need to be exported later
                event.source_frag = "Removed"


class Event():
    '''represents the symbols, also known as the "events," within a SV transformation'''

    def __init__(self, sv_parent, length, length_range, symbol, source_frag=None):
        '''
        sv_parent: Structural Variant, event is always part of larger SV
        '''
        self.sv_parent = sv_parent
        self.length = length
        self.length_range = length_range
        self.symbol = symbol  # refers to symbol in SV's transformation
        self.source_chr = None if not source_frag else source_frag
        self.source_frag = None
        self.start = None
        self.end = None

    def __repr__(self):
        return "<Event {}>".format({"length": self.length, "symbol": self.symbol, "start": self.start, "end": self.end,
                                    "source_chr": self.source_chr, "source_frag": self.source_frag})


class Blocks():
    """
    Groups together source and target symbols between dispersion events (_)
    """
    def __init__(self, sv):
        self.sv = sv
        self.source_blocks = []
        self.target_blocks = []
        self.generate_blocks()
        # TODO: optional dispersion flip should be done in init step
        # if _____:
        #   self.dispersion_flip()
        self.track_original_symbol()

    def generate_blocks(self):
        # *** DO WE EVER USE SOURCE BLOCKS?? DON'T WE JUST NEED THE TARGET BLOCKS WITH THE RIGHT REF FRAGS?
        self.source_blocks = self.find_blocks(self.sv.source_unique_char)
        self.target_blocks = self.find_blocks(self.sv.target_unique_char)

    def find_blocks(self, transformation):
        # transformation: tuple of strings (source or target chars)
        # -> returns list of lists of subevents
        # Ex. ("A","B","_","C","D") -> [[Event("A",...),Event("B",...)], [Event("_1")], [Event("C",...),Event("D",...)]]
        # Ex. ("A","B","_","_") -> [[Event("A",...),Event("B",...)],[Event("_1")],[Event("_2")]]
        blocks = [[]]
        for symbol in transformation:
            if symbol.startswith(Symbols.DIS.value):
                # going to add singleton lists with the dispersions where they occur so we can keep track of the sizes
                source_event = self.sv.events_dict[symbol]
                disp_event = Event(sv_parent=self.sv, length=source_event.length, length_range=None, symbol=symbol)
                blocks.append([disp_event])
                blocks.append([])
            else:
                # used to find corresponding event from encoding, all keys in encoding are in uppercase
                source_event = self.sv.events_dict[symbol[0].upper()]
                target_event = Event(sv_parent=self.sv, length=source_event.length, length_range=None, symbol=symbol)
                blocks[-1].append(target_event)
        return blocks

    @staticmethod

    def dispersion_flip(self):
        # perform the optional flip of the blocks list dictated by the flipped dispersion
        # --> caveats: - how to determine which section of the source/target blocks lists to flip?
        # -->          - how to trigger this/identify or mark a dispersion to be flipped?
        # -->          - should it be a quality of a whole dispersion-related event?
        pass

    def track_original_symbol(self):
        # finds which region/block the original symbol is in
        # if symbol is later found in another region, then translocation detected
        # blocks: list of lists of symbols
        for idx, block in enumerate(self.source_blocks):
            for event in block:
                if len(event.symbol) == 1:  # means it's an original symbol
                    event.original_block_idx = idx

    def __repr__(self):
        return f"SOURCE_BLOCKS: {self.source_blocks}\n" \
               f"TARGET_BLOCKS: {self.target_blocks}"
