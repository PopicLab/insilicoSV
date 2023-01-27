from constants import *
import utils
import random


class Structural_Variant():
    def __init__(self, sv_type, mode, length_ranges=None, source=None, target=None, vcf_rec=None, ref_fasta=None):
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
        vcf_rec: (to be used in fixed mode) vcf record giving sv information that will instantiate the event
        ref_fasta: for extracting reference frag for vcf records in fixed mode initialization
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
        # *** are these sv.start/end attributes even used?
        self.start = None  # defines the space in which SV operates
        self.end = None
        self.req_space = None  # required space for SV, sum of event lengths
        # TODO: only used for iterating over events in choose_rand_pos() -- delete and replace with source_blocks?
        self.source_events = []  # list of Event classes for every symbol in source sequence
        self.events_dict = dict()  # maps every unique symbol in source and target to an Event class
        self.changed_fragments = []  # list recording new fragments to be placed in the output ref
        if self.type in [Variant_Type.dDUP, Variant_Type.INV_dDUP, Variant_Type.div_dDUP, Variant_Type.TRA]:
            # if not self.dispersion_flip and random.randint(0, 1):
            if True:
                self.dispersion_flip = True
        # initialize_events sets the values of events_dict, source_dict, and req_space
        if mode == 'randomized':
            self.initialize_events(length_ranges)
        else:
            # ** in the case of fixed mode, the sim.process_vcf() method seems to act as the alternative initialization
            # ** method -- better practice to have a fixed mode initialization method here that we'll feed the vcf record to?
            self.initialize_events_fixed(vcf_rec, ref_fasta)

        sv_blocks = Blocks(self)
        # debug
        # print(f'==== sv_blocks object ====\n{sv_blocks}')
        self.source_symbol_blocks = sv_blocks.source_blocks
        self.target_symbol_blocks = sv_blocks.target_blocks

        if mode == 'randomized':
            # specifies if sv is unable to be simulated due to random placement issues
            # will be turned on later
            self.active = False

            # 1 = homozygous, 0 = heterozygous
            self.ishomozygous = Zygosity.UNDEFINED

            # stores list of booleans specifying if SV will be applied to certain haplotype for assigned chromosome
            self.hap = [False, False]
        else:
            # if in fixed mode, will need to manually call assign_locations (since that is called by choose_rand_pos() in random mode)
            self.assign_locations(self.start)

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

        # TODO: need to fix the way events in events_dict have their source_frags assigned (problem: the source
        #  frag is being assigned with respect to the forward-facing subevent orientation whether or not our
        #  SV is a flipped-dispersion event -- then in the case of a flipped-dispersion event, the source frags
        #  taken from the dict wont be accurate)
        #    ---> is it correct to just flip the order of source_events if our sv is meant to be flipped? (by the
        #         source_frag assignment in choose_rand_pos, it looks like it's jsut determined by order in the list)
        if self.dispersion_flip:
            self.source_events = self.source_events[::-1]

        self.req_space = sum([event.length for event in self.source_events])

        # # debug
        # print('END OF INITIALIZE_EVENTS')
        # print('sv.events_dict:')
        # for ev in self.events_dict.keys():
        #     print(self.events_dict[ev])

    def initialize_events_fixed(self, vcf_record, ref_fasta):
        """
        initialization method for SV read in from vcf -- need to prepare objects needed for creation of target blocks
        --> needs to fully populate events_dict with start/end (reflecting optional dispersion flip), lengths, and source frags
        --> (with those in place, Blocks() and assign_locations() should have everything they need to run before change_frags())
        """
        # debug
        print('==== INITIALIZE_EVENTS (FIXED-MODE) ====')
        print(f'source symbols : {self.source_unique_char}')
        print(f'target symbols : {self.target_unique_char}')
        # TODO: populate events_dict with information from vcf record
        #  --> what's the best way to handle events with multiple input symbols?
        #  ----> anything like that will have multiple vcf records, but not sure how to write that logic here when
        #  ----> we're just looking at a single record
        source_len = vcf_record.stop - vcf_record.start if 'SVLEN' not in vcf_record.info else vcf_record.info['SVLEN']
        for symbol in self.source_unique_char:
            if symbol == 'A':
                source_ev = Event(self, symbol, source_len, (source_len, source_len), symbol)
                source_ev.start = vcf_record.start
                source_ev.end = vcf_record.stop
                source_ev.source_chr = vcf_record.chrom
                source_ev.source_frag = ref_fasta.fetch(source_ev.source_chr, source_ev.start, source_ev.end)
                self.events_dict[symbol] = source_ev
            if symbol.startswith(Symbols.DIS.value):
                flipped_dispersion = vcf_record.info['TARGET'] < vcf_record.start
                disp_len = vcf_record.info['TARGET'] - vcf_record.stop if not flipped_dispersion else \
                    vcf_record.start - vcf_record.info['TARGET']
                disp_ev = Event(self, symbol, disp_len, (disp_len, disp_len), symbol)
                disp_ev.start = vcf_record.stop if not flipped_dispersion else vcf_record.info['TARGET']
                disp_ev.end = vcf_record.info['TARGET'] if not flipped_dispersion else vcf_record.start
                disp_ev.source_chr = vcf_record.chrom
                disp_ev.source_frag = ref_fasta.fetch(disp_ev.source_chr, disp_ev.start, disp_ev.end)
                self.events_dict[symbol] = disp_ev
        # debug
        print('events_dict = ')
        for ev in self.events_dict:
            print(ev)

        # TODO: populate other SV attributes (everything populated in current process_vcf logic)
        pass

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
        # # debug
        # print('===LOCATIONS NOT ASSIGNED YET===\ntarget_symbol_blocks:')
        # for bl in self.target_symbol_blocks:
        #     print(bl)
        # print('events_dict:')
        # for ev in self.events_dict.keys():
        #     print(self.events_dict[ev])

        # Trying logic based on the position assignment of source events in choose_rand_pos()
        # (now that we're flipping the sv.source_events list for flipped-dispersion events)
        for block in self.target_symbol_blocks:
            for ev in block:
                # if the event is one also found in the source, place it at the location given in events_dict
                # --> the events that stay the same will need to be in the same place in both input and output ref
                if ev.symbol.upper() in self.events_dict.keys():
                    source_ev = self.events_dict[ev.symbol.upper()]
                    ev.start = source_ev.start
                    ev.end = source_ev.end
                    ev.source_frag = self.get_event_frag(source_ev, ev.symbol)
        # ===> then everything that wasn't just assigned will be a new event (all to be modeled as insertion fragments)
        # ------> thus will just need to have start/end set to the nearest event boundary from the ones placed above^
        # --> position is just determined by event adjacency, easier to ignore block boundaries here
        target_events = [ev for bl in self.target_symbol_blocks for ev in bl]
        # singleton event that's novel (i.e., INS)
        if len(target_events) == 1 and target_events[0].start is None:
            ev = target_events[0]
            ev.start = start_pos
            ev.end = start_pos
            # position assigned, need to get source frag
            source_event = self.events_dict[ev.symbol[0].upper()]
            ev.source_frag = self.get_event_frag(source_event, ev.symbol)
        else:
            for i in range(len(target_events)):
                ev = target_events[i]
                if ev.start is None:
                    if i == 0:
                        # if the first event is novel, set start/end to the start of the nearest event
                        j = i + 1
                        while target_events[j].start is None:
                            j += 1
                        ev.start = target_events[j].start
                        ev.end = target_events[j].start
                    else:
                        ev.start = target_events[i - 1].end
                        ev.end = target_events[i - 1].end
                    # position assigned, need to get source frag
                    source_event = self.events_dict[ev.symbol[0].upper()]
                    ev.source_frag = self.get_event_frag(source_event, ev.symbol)

        # debug
        # print('===LOCATIONS ASSIGNED===\ntarget_symbol_blocks:')
        # for bl in self.target_symbol_blocks:
        #     print(bl)

    def change_fragment(self):
        '''
        Takes the mapping of symbols to events and the target sequence to construct a replacement sequence for the reference fragment
        '''
        changed_fragments = []
        assert (self.start is not None and self.end is not None), "Undefined SV start for {}".format(
            self)  # start & end should have been defined alongside event positions
        block_start = None
        block_end = None

        # # debug
        print('===CHANGE_FRAGMENT===')
        print('sv.events_dict:')
        for ev in self.events_dict.keys():
            print(self.events_dict[ev])
        print(f'target blocks = {self.target_symbol_blocks}')
        # special case: deletion -- len(target_symbol_blocks) == 0
        if self.target_symbol_blocks == [[]]:
            changed_fragments.append([self.start_chr, self.start, self.end, ''])
        else:
            for idx, block in enumerate(self.target_symbol_blocks):
                new_frag = ''
                if len(block) == 0:
                    # this branch will be executed for TRAs:
                    # --> want to delete the A-length interval on the opposite side of the dispersion as our A'
                    if idx == 0:
                        del_len = len(self.target_symbol_blocks[2][0].source_frag)
                        disp_ev = self.target_symbol_blocks[1][0]
                        block_start = disp_ev.start - del_len
                        block_end = disp_ev.start
                    else:
                        del_len = len(self.target_symbol_blocks[idx - 2][0].source_frag)
                        disp_ev = self.target_symbol_blocks[idx - 1][0]
                        block_start = disp_ev.end
                        block_end = block_start + del_len
                    changed_fragments.append([self.start_chr, block_start, block_end, new_frag])
                    continue
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
                # print(f'new change fragment : {changed_fragments[-1]}')

        self.changed_fragments = changed_fragments
        # debug
        print(f'=== CHANGED_FRAGMENTS FOR {self.type} ===\n{self.changed_fragments}')
        self.clean_event_storage()  # clean up unused storage - we do not need to store most source_frags anymore
        return changed_fragments

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
        self.source_chr = None
        self.source_frag = None if not source_frag else source_frag
        self.start = None
        self.end = None

    def __repr__(self):
        return "<Event {}>".format({"length": self.length, "symbol": self.symbol, "start": self.start, "end": self.end,
                                    "source_chr": self.source_chr,
                                    "source_frag": self.source_frag if not self.symbol.startswith(Symbols.DIS.value) else
                                    'frag omitted'})


class Blocks():
    """
    Groups together source and target symbols between dispersion events (_)
    """
    def __init__(self, sv):
        self.sv = sv
        self.source_blocks = []
        self.target_blocks = []
        self.generate_blocks()
        # optional dispersion flip should be done in init step
        if self.sv.type in [Variant_Type.TRA, Variant_Type.dDUP, Variant_Type.INV_dDUP, Variant_Type.div_dDUP] \
                and self.sv.dispersion_flip:
            self.flip_blocks()
        self.track_original_symbol()

    def generate_blocks(self):
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
                # if symbol different from source symbol then being added to input ref
                if symbol.upper() != source_event.symbol:
                    target_event.length = 0
                blocks[-1].append(target_event)
        return blocks

    def flip_blocks(self):
        # perform the optional flip of the blocks list dictated by the flipped dispersion
        # --> caveats: - how to determine which section of the source/target blocks lists to flip?
        # -->          - how to trigger this/identify or mark a dispersion to be flipped?
        # -->          - should it be a quality of a whole dispersion-related event?
        # ===> Is it incorrect to just flip the entire list of target blocks? I don't think this is the completely
        # general rule, but for all of our dispersion events I think this will always be the case
        self.source_blocks = self.source_blocks[::-1]
        self.target_blocks = self.target_blocks[::-1]

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
