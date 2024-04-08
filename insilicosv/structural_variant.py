import random

from insilicosv import utils
from insilicosv.constants import \
    SV_KEY, Symbols, Variant_Type, Zygosity, DISPERSION_TYPES

class Structural_Variant:
    def __init__(self, sv_type, mode, length_ranges=None, source=None, target=None, vcf_rec=None, ref_fasta=None,
                 overlap_event=None, div_prob=None):
        """
        sv_type: Enum either specifying one of the prewritten classes or a Custom transformation
        mode: flag to indicate whether constructor called in fixed or randomized mode
        length_ranges: list containing tuple(s) (min_length, max_length) or singleton int
        source: tuple representing source sequence, optional
        target: tuple representing target sequence, optional
        vcf_rec: (fixed mode) vcf record giving sv information that will instantiate the event
        ref_fasta: for extracting reference frag for vcf records in fixed mode initialization
        overlap_event: (chr, start, end, elt_type) tuple representing the genome element that this SV is meant to overlap, optional
        """
        self.type = sv_type
        if self.type != Variant_Type.Custom:
            self.source, self.target = SV_KEY[self.type]
            self.name = self.type.name
        else:
            self.source, self.target = tuple(source), tuple(target)
            self.name = str("".join(self.source)) + ">" + str("".join(self.target))

        utils.validate_symbols(self.source, self.target)
        self.source_unique_char, self.target_unique_char = \
            Structural_Variant.add_unique_dispersion_ids(self.source), \
            Structural_Variant.add_unique_dispersion_ids(self.target)

        self.start = None
        self.end = None
        self.start_chr = None
        self.req_space = None  # required space for SV, sum of event lengths
        self.source_events = []  # list of Event classes for every symbol in source sequence
        self.events_dict = dict()  # maps every unique symbol in source and target to an Event class
        self.changed_fragments = []  # list recording new fragments to be placed in the output ref
        self.dispersion_flip = False  # orientation of dispersion event
        self.insseq_from_rec = None  # space to store INSSEQ for fixed-mode INS event
        self.overlap_event = overlap_event  # <- element tuple of form (chrom, start, end, type) (optional)
        self.div_prob = div_prob  # <- divergence probability param (optionally given for type == DIVERGENCE)

        if self.type in DISPERSION_TYPES:
            if random.randint(0, 1):
                self.dispersion_flip = True
        if mode == 'randomized':
            self.initialize_events(length_ranges)
        else:
            self.initialize_events_fixed(vcf_rec, ref_fasta)

        self.sv_blocks = Blocks(self)
        self.target_symbol_blocks = self.sv_blocks.target_blocks

        if mode == 'randomized':
            self.active = False
            self.ishomozygous = Zygosity.UNDEFINED
            self.hap = [False, False]
        else:
            # manual call to assign_locations in fixed mode
            self.assign_locations(self.start)

    def __repr__(self):
        return "<SV transformation \"{}\" -> \"{}\" taking up {} non-dispersion spaces>".format(
            ''.join(self.source), ''.join(self.target),
            sum([event.length for event in self.source_events if not event.symbol.startswith(Symbols.DIS.value)]))

    @staticmethod
    def add_unique_dispersion_ids(transformation):
        """In a tuple defining the source or target of a transformation (see constants.SV_KEY),
        make each dispersion event unique by appending a unique ID.

        Examples:
           ("A", "_") becomes ("A", "_1")
           ("A", "_", "B", "_", "C") becomes ("A", "_1", "B", "_2", "C")
        """
        unique_transform = []
        unique_id = 1
        for component in transformation:
            if component != Symbols.DIS.value:
                unique_transform.append(component)
            else:
                unique_transform.append(component + str(unique_id))
                unique_id += 1
        return tuple(unique_transform)

    def get_event_frag(self, event, symbol):
        """Returns the target event sequence obtained from the source event
        sequence by applying the specified target transformation (if specified).
        
        For example, if the symbol specifies that the source sequence
        should be complemented or mutated, that transformation is applied.
        
        Args:
          event: source event from events_dict
          symbol: target symbol
        
        Returns:
          The transformed sequence.
        """
        assert 1 <= len(symbol) <= 2
        assert symbol[0].isalpha() or symbol[0] == Symbols.DIS.value

        decode_funcs = {
            "invert": lambda string: utils.complement(string[::-1]),
            "diverge": lambda string: utils.divergence(string,
                                                       divergence_prob=(1 if self.type == Variant_Type.SNP
                                                                        else self.div_prob))
        }

        if symbol[0].islower():
            return decode_funcs["invert"](event.source_frag)
        elif symbol[-1] == Symbols.DIV.value:
            return decode_funcs["diverge"](event.source_frag)
        else:
            return event.source_frag

    def initialize_events(self, length_ranges):
        """
        Initializes event classes and creates a mapping of symbol to event

        length_ranges: list of tuples specifying min and max length for events within SV
        -> returns list of events in source sequence
        """
        all_symbols = []
        for ele in self.source_unique_char + self.target_unique_char:
            if len(ele) > 0 and (len(ele) == 1 or ele.startswith(Symbols.DIS.value)) and ele.upper() not in all_symbols:
                all_symbols.append(ele.upper())
        all_symbols.sort()

        # symbol_len: (key = symbol, value = chosen length)
        symbol_len = dict()
        if len(length_ranges) > 1:  # values given by user represents custom ranges for each event symbol of variant in lexicographical order
            assert (len(length_ranges) == len(all_symbols)), \
                "Number of length_ranges entered does not match the number of symbols (remember foreign insertions and dispersions) present!"
            for idx, symbol in enumerate(all_symbols):
                symbol_len[symbol] = random.randint(length_ranges[idx][0], length_ranges[idx][1])
        elif len(length_ranges) == 1:  # value given by user represents length (same range) of each event within variant in lexicographical order
            for symbol in all_symbols:
                symbol_len[symbol] = random.randint(length_ranges[0][0], length_ranges[0][1])
        else:
            raise Exception("length_ranges parameter expects at least one tuple")

        ovlp_frag = random.choice([frag for frag in self.source_unique_char if frag[0] != '_']) if self.overlap_event is not None else None
        for idx, symbol in enumerate(all_symbols):
            if self.overlap_event is not None and symbol == ovlp_frag:
                ovlp_event_len = int(self.overlap_event[2]) - int(self.overlap_event[1])
                event = Event(length=ovlp_event_len, symbol=symbol,
                              start=int(self.overlap_event[1]), end=int(self.overlap_event[2]))
            else:
                event = Event(length=symbol_len[symbol], symbol=symbol)
            self.events_dict[symbol] = event

        for symbol in self.source_unique_char:
            self.source_events.append(self.events_dict[symbol])

        if self.dispersion_flip:
            self.source_events = self.source_events[::-1]

        self.req_space = sum([event.length for event in self.source_events])

    def initialize_events_fixed(self, vcf_record, ref_fasta):
        # initialization method for SV read in from vcf
        source_len = vcf_record.stop - vcf_record.start if 'SVLEN' not in vcf_record.info else vcf_record.info['SVLEN']
        for symbol in self.source_unique_char:
            if symbol == Symbols.REQUIRED_SOURCE.value:
                source_ev = Event(length=source_len, symbol=symbol)
                source_ev.start = vcf_record.start
                source_ev.end = vcf_record.stop
                source_ev.source_chr = vcf_record.chrom
                source_ev.source_frag = ref_fasta.fetch(source_ev.source_chr, source_ev.start, source_ev.end) if \
                    vcf_record.id != Variant_Type.SNP.value else vcf_record.ref
                self.events_dict[symbol] = source_ev
            if symbol.startswith(Symbols.DIS.value):
                self.dispersion_flip = vcf_record.info['TARGET'] < vcf_record.start
                disp_len = vcf_record.info['TARGET'] - vcf_record.stop if not self.dispersion_flip else \
                    vcf_record.start - vcf_record.info['TARGET']
                disp_ev = Event(length=disp_len, symbol=symbol)
                disp_ev.start = vcf_record.stop if not self.dispersion_flip else vcf_record.info['TARGET']
                disp_ev.end = vcf_record.info['TARGET'] if not self.dispersion_flip else vcf_record.start
                disp_ev.source_chr = vcf_record.chrom
                disp_ev.source_frag = ref_fasta.fetch(disp_ev.source_chr, disp_ev.start, disp_ev.end)
                self.events_dict[symbol] = disp_ev
        # storing the insertion seq for INSs/SNPs with an insertion sequence / alt allele given in the vcf
        if vcf_record.id == 'SNP':
            # NB: only supporting SNP records with a single allele ALT reported
            self.insseq_from_rec = vcf_record.alts[0]
        if 'SVTYPE' in vcf_record.info and vcf_record.info['SVTYPE'] == 'INS':
            if 'INSSEQ' in vcf_record.info:
                self.insseq_from_rec = vcf_record.info['INSSEQ'][0]
            source_ev = Event(length=source_len, symbol=Symbols.REQUIRED_SOURCE.value)
            self.events_dict[Symbols.REQUIRED_SOURCE.value] = source_ev
            self.start = vcf_record.start
            self.end = self.start
            print(f'======= TEST {self.start=}, {self.end=}')
        else:
            self.start = self.events_dict[Symbols.REQUIRED_SOURCE.value].start
            self.end = self.events_dict[Symbols.REQUIRED_SOURCE.value].end if '_1' not in self.events_dict else self.events_dict['_1'].end
        self.start_chr = vcf_record.chrom

        # handling for divergent repeat simulation logic (div_dDUPs placed into R1 need to correspond to dDUPs in R2)
        if self.type == Variant_Type.div_dDUP:
            self.target_unique_char = ("A", "_1", "A'")

        self.active = True
        if vcf_record.samples[0]['GT'] == (1, 1):
            self.ishomozygous = Zygosity.HOMOZYGOUS
            self.hap = [True, True]
        else:
            self.ishomozygous = Zygosity.HETEROZYGOUS
            self.hap = random.choice([[True, False], [False, True]])

    def assign_locations(self, start_pos):
        """
        assign events start and end positions (once target blocks are populated and in the right order)
        """
        for block in self.target_symbol_blocks:
            for ev in block:
                ev.source_chr = self.start_chr
                # if the event is one also found in the source, place it at the location given in events_dict
                if ev.symbol.upper() in self.events_dict or ev.symbol[-1] == Symbols.DIV.value:
                    source_ev = self.events_dict[ev.symbol.upper() if ev.symbol.upper() in self.events_dict else ev.symbol[0]]
                    ev.start = source_ev.start
                    ev.end = source_ev.end
                    ev.source_frag = self.get_event_frag(source_ev, ev.symbol) if self.insseq_from_rec is None else self.insseq_from_rec

        # everything that wasn't assigned above will be modeled as insertion fragment placed at the nearest event boundary
        target_events = [ev for bl in self.target_symbol_blocks for ev in bl]
        # singleton event that's novel (i.e., INS)
        if len(target_events) == 1 and target_events[0].start is None:
            ev = target_events[0]
            ev.start = start_pos
            ev.end = start_pos
            source_event = self.events_dict[ev.symbol[0].upper()]
            if self.insseq_from_rec is not None:
                ev.source_frag = self.insseq_from_rec
            else:
                ev.source_frag = self.get_event_frag(source_event, ev.symbol)
                if ev.source_frag is None:
                    ev.source_frag = utils.generate_seq(ev.length)
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
                    source_event = self.events_dict[ev.symbol[0].upper()]
                    ev.source_frag = self.get_event_frag(source_event, ev.symbol)

        self.sv_blocks.generate_target_events_dict()

    def change_fragment(self):
        """
        Takes the mapping of symbols to events and the target sequence to construct a replacement sequence for the reference fragment
        """
        changed_fragments = []
        assert (self.start is not None and self.end is not None), "Undefined SV start for {}".format(
            self)
        block_start = None
        block_end = None

        if self.target_symbol_blocks == [[]]:  # special case: simple deletion -- len(target_symbol_blocks) == 0
            changed_fragments.append([self.start_chr, self.start, self.end, ''])
        else:
            for idx, block in enumerate(self.target_symbol_blocks):
                new_frag = ''
                if len(block) == 0 or block[0].symbol.startswith(Symbols.DIS.value):
                    continue
                for i in range(len(block)):
                    ev = block[i]
                    new_frag += ev.source_frag
                    if i == 0:
                        block_start = ev.start
                    if i == len(block) - 1:
                        block_end = ev.end
                changed_fragments.append([self.start_chr, block_start, block_end, new_frag])
            # Create a deletion over the source symbol unless it appears in the target as an INV or DIV
            target_symbols = [ev.symbol.upper() for bl in self.target_symbol_blocks for ev in bl]
            for source_sym in self.events_dict.keys():
                if source_sym not in target_symbols and source_sym + Symbols.DIV.value not in target_symbols:
                    del_ev = self.events_dict[source_sym]
                    changed_fragments.append([del_ev.source_chr, del_ev.start, del_ev.end, ''])

        self.changed_fragments = changed_fragments
        self.clean_event_storage()

    def clean_event_storage(self):
        # remove source fragments from events to save space as they are no longer needed
        for event in self.events_dict.values():
            if event.symbol in self.source_unique_char:  # do not clean out insertion fragments as they'll need to be exported later
                event.source_frag = "Removed"


class Event:
    """
    represents the symbols, also known as the "events," within a SV transformation
    """

    def __init__(self, *, length, symbol, source_frag=None, start=None, end=None):
        self.length = length
        self.symbol = symbol  # refers to symbol in SV's transformation
        self.source_chr = None
        self.source_frag = None if not source_frag else source_frag
        self.start = start
        self.end = end

    def __repr__(self):
        return "<Event {}>".format({"length": self.length, "symbol": self.symbol, "start": self.start, "end": self.end,
                                    "source_chr": self.source_chr})


class Blocks:
    """
    Groups together target symbols between dispersion events
    """
    def __init__(self, sv):
        self.sv = sv
        self.target_blocks = []
        self.generate_blocks()
        if self.sv.type in DISPERSION_TYPES and self.sv.dispersion_flip:
            self.flip_blocks()
        self.target_events_dict = None

    def generate_blocks(self):
        self.target_blocks = self.find_blocks(self.sv.target_unique_char)

    def find_blocks(self, transformation):
        # transformation: tuple of strings (source or target chars)
        # -> returns list of lists of events
        # Ex. ("A","B","_","C","D") -> [[Event(symbol="A",...),Event(symbol="B",...)],
        #                               [Event(symbol="_1")],
        #                               [Event(symbol="C",...),Event(symbol="D",...)]]
        # Ex. ("A","B","_","_") -> [[Event(symbol="A",...),Event(symbol="B",...)],
        #                           [Event(symbol="_1")],
        #                           [Event(symbol="_2")]]
        blocks = [[]]
        for symbol in transformation:
            if symbol.startswith(Symbols.DIS.value):
                source_event = self.sv.events_dict[symbol]
                disp_event = Event(length=source_event.length, symbol=symbol)
                blocks.append([disp_event])
                blocks.append([])
            else:
                # used to find corresponding event from encoding, all keys in encoding are in uppercase
                source_event = self.sv.events_dict[symbol[0].upper()]
                target_event = Event(length=source_event.length, symbol=symbol)
                # if symbol different from source symbol then being added to input ref
                if symbol.upper() != source_event.symbol:
                    target_event.length = 0
                blocks[-1].append(target_event)
        return blocks

    def flip_blocks(self):
        # perform the optional flip of the blocks list dictated by the flipped dispersion
        self.target_blocks = self.target_blocks[::-1]

    def generate_target_events_dict(self):
        # setter for target_events_dict attribute (to be populated after location assignment of target events)
        self.target_events_dict = {ev.symbol: ev for b in self.target_blocks for ev in b}

    def __repr__(self):
        return f"TARGET_BLOCKS: {self.target_blocks}"

