from collections import defaultdict
from contextlib import closing, contextmanager
import dataclasses
from dataclasses import dataclass, field
from enum import Enum
import functools
import itertools
import logging
import operator
import os
import os.path
import random
import re
from typing_extensions import TypeAlias, Optional, Any, Union, Iterable

from intervaltree import IntervalTree, Interval
import pysam
from sortedcontainers import SortedSet  # type: ignore

logger = logging.getLogger(__name__)

def if_not_none(a, b):
    return a if a is not None else b

def chk(cond, msg="Error") -> None:
    if not cond:
        raise RuntimeError(f"insilicoSV error: {msg}")

def generate_seq(length: int) -> str:
    base_map = {1: "A", 2: "T", 3: "G", 4: "C"}
    return ''.join([base_map[random.randint(1, 4)] for _ in range(length)])

def complement(seq: str) -> str:
    output = ""
    base_complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    for base in seq.upper():
        if base in base_complements:
            output += base_complements[base]
        else:
            output += base

    return output

def reverse_complement(seq: str) -> str:
    return complement(seq[::-1])

def remove_file(path: str) -> None:
    if os.path.exists(path):
        os.remove(path)

def divergence(seq: str, divergence_prob: float) -> str:
    # apply random base flips to input sequence
    p = float(divergence_prob)
    assert 0 <= p <= 1
    if p == 0:
        return seq
    def mutate_base(b: str) -> str:
        if b not in 'TCGA':
            return b
        return random.choice(list({"A", "C", "T", "G"} - {b}))
    return ''.join([b if random.random() > p else mutate_base(b) for b in seq.upper()])

def is_readable_file(fname: str) -> bool:
    try:
        with open(fname) as test_f:
            pass
        return True
    except:
        return False

def as_list(val: Union[str, list[str]]) -> list[str]:
    if isinstance(val, str):
        val = [val]
    return val

def is_valid_closed_int_range(val) -> bool:
    return (
        isinstance(val, list) and
        len(val) == 2 and
        isinstance(val[0], int) and
        isinstance(val[1], int) and
        0 <= val[0] <= val[1])

def is_list_of(type_, val) -> bool:
    return isinstance(val, list) and all(isinstance(v, type_) for v in val)

def is_convertible_to(type_, val) -> bool:
    try:
        _ = type_(val)
        return True
    except Exception:
        return False

@dataclass(order=True, frozen=True)
class Locus:
    chrom: str
    pos: int

    def shifted(self, delta: int) -> 'Locus':
        return Locus(chrom=self.chrom, pos=self.pos + delta)

    def shifted_to(self, pos: int) -> 'Locus':
        return Locus(chrom=self.chrom, pos=pos)

@dataclass(order=True, frozen=True)
class Region:
    """A specific, contiguous genomic region: (chrom, start, end).
    
    The region includes `start` but not `end`.  An empty region with start==end represents
    an insertion point before `start`.  Insertions at the same point are ordered by `order_key`.
    The coordinates `start` and `end` are 0-based.
    """

    chrom: str
    start: int
    end: int
    order_key: tuple = ()

    # Optional extra info that can be associated with Regions representing
    # ROIs or pieces of ROIs:

    # region type (e.g. L1, Alu, etc)
    kind: str = ''

    # additional per-region data, such as repeat unit info for tandem repeats
    data: str = ''

    # if this region is derived from an ROI, start/end of the original ROI
    orig_start: int = -1
    orig_end: int = -1

    def __post_init__(self) -> None:
        assert self.chrom
        assert self.start <= self.end

    def is_empty(self) -> bool:
        return self.start == self.end

    def length(self) -> int:
        return self.end - self.start

    def replace(self, **kw) -> 'Region':
        """Return a copy of self with fields replaced according to `kw`"""
        return dataclasses.replace(self, **kw)

    def head(self, length: int) -> 'Region':
        assert 0 <= length <= self.length()
        return self.replace(end=self.start + length)
    
    def tail(self, length: int) -> 'Region':
        assert 0 <= length <= self.length()
        return self.replace(start=self.end - length)

    def shifted(self, delta: int) -> 'Region':
        return self.replace(start=self.start + delta, end=self.end + delta)

    def padded(self, padding: int) -> 'Region':
        return self.replace(start=self.start - padding, end=self.end + padding)

    def shifted_to(self, *, start=None, end=None) -> 'Region':
        assert (start is None) != (end is None)
        if start is not None:
            return self.replace(start=start, end=start+self.length())
        else:
            return self.replace(start=end-self.length(), end=end)

    def is_adjacent_to(self, other: 'Region') -> bool:
        return (self.chrom == other.chrom and
                ((self.end == other.start) or
                 (self.start == other.end)))

# end: class Region

@dataclass(frozen=True)
class RegionFilter:
    region_kinds: Optional[tuple[str, ...]] = None
    region_length_range: tuple[Optional[int], Optional[int]] = (None, None)

    def satisfied_for(self, region: Region) -> bool:
        if (self.region_kinds is not None and
            (not region.kind or
             not any((region_kinds.upper() == 'ALL' and region.kind != '_reference_') or
                     region.kind.startswith(region_kinds)
                     for region_kinds in self.region_kinds))):
            return False
        if (self.region_length_range[0] is not None and 
            region.length() < self.region_length_range[0]):
            return False
        if (self.region_length_range[1] is not None and
            region.length() > self.region_length_range[1]):
            return False
        return True

class OverlapMode(Enum):
    
    """Types of spatial relationship between an anchor (SV sub-region) and an ROI."""

    PARTIAL = "partial"      # exactly one endpoint of anchor in ROI
    CONTAINED = "contained"  # both endpoints of anchor in ROI
    EXACT = "exact"  # anchor exactly matches the ROI

def n_valid_placements(*, overlap_mode: OverlapMode, region_length: int, anchor_length: int) -> int: 
    if overlap_mode == OverlapMode.CONTAINED:
        return max(0, region_length - anchor_length + (1 if anchor_length > 0 else -1))
    if overlap_mode == OverlapMode.PARTIAL:
        return 2 * min(max(0, anchor_length - 1), max(0, region_length - 1))
    assert False

class RegionSet:
    """A collection of genomic regions"""

    chrom2itree: dict[str, IntervalTree]
    # IntervalTree does not support empty intervals,
    # so we store them separately
    chrom2empty_regions: dict[str, SortedSet]  

    def __init__(self, regions=None):
        regions = regions or []
        chrom2nonempty_regions = defaultdict(list)
        self.chrom2empty_regions = defaultdict(
            lambda: SortedSet(key=lambda region: region.start))
        for region in regions:
            if region.is_empty():
                self.chrom2empty_regions[region.chrom].add(region)
            else:
                chrom2nonempty_regions[region.chrom].append(region)

        self.chrom2itree = defaultdict(IntervalTree)
        for chrom, chrom_nonempty_regions in chrom2nonempty_regions.items():
            if len(chrom_nonempty_regions) > 100000:
                logger.debug(f'RegionSet init: {chrom=} {len(chrom_nonempty_regions)=}')
            self.chrom2itree[chrom] = IntervalTree.from_tuples((region.start, region.end, region)
                                                         for region in chrom_nonempty_regions)

    def assert_valid(self) -> None:
        for chrom, chrom_itree in self.chrom2itree.items():
            for ival in chrom_itree:
                assert isinstance(ival.data, Region)
                assert chrom == ival.data.chrom
                assert ival.begin == ival.data.start
                assert ival.end == ival.data.end
                assert not ival.data.is_empty()

        for chrom, chrom_empty_regions in self.chrom2empty_regions.items():
            for region in chrom_empty_regions:
                assert region.chrom == chrom
                assert region.is_empty()
            assert all(r1 <= r2 for r1, r2 in pairwise(chrom_empty_regions))

    def __len__(self) -> int:
        return (sum(len(chrom_itree) for chrom_itree in self.chrom2itree.values()) +
                sum(len(chrom_empty_regions) for chrom_empty_regions in self.chrom2empty_regions.values()))

    def __contains__(self, region: Region) -> bool:
        return ((not region.is_empty() and 
                 Interval(region.start, region.end, region) in self.chrom2itree[region.chrom]) or
                (region.is_empty() and region in self.chrom2empty_regions[region.chrom]))
                
    @staticmethod
    def from_beds(bed_paths: list[str], verbose: bool = False) -> 'RegionSet':
        regions = []
        for bed_path in bed_paths:
            logger.info(f'reading bed file {bed_path}')
            with open(bed_path) as bed:
                for line_num, line in enumerate(bed):
                    if (verbose and (line_num % 500000) == 0):
                        logger.debug(f'line {line_num}')
                    loc = f'file {bed_path}, line {line_num}'
                    if line.startswith('#') or line.isspace():
                        continue
                    fields = line.strip().split()
                    chk(len(fields) >= 4,
                        f'{loc}: too fiew fields in line')
                    chrom, start_str, end_str, kind = fields[:4]
                    chk(all((chrom, start_str, end_str, kind)),
                        f'{loc}: empty value in first four columns')
                    try:
                        start, end = int(start_str), int(end_str)
                    except ValueError:
                        chk(False, f'{loc}: invalid start or end of region')
                    chk(start < end, f'{loc}: region start must be less than end')
                    chk(0 <= start, f'{loc}: region start must be non-negative')

                    data = fields[4] if len(fields) >= 5 else ''

                    regions.append(Region(chrom=chrom, start=start,
                                          end=end, kind=kind, data=data,
                                          orig_start=start, orig_end=end))

        logger.debug(f'Constructing RegionSet from {len(regions)} regions...')
        region_set = RegionSet(regions)
        logger.debug(f'Constructed RegionSet from {len(regions)} regions.')
        return region_set

    @staticmethod
    def from_vcf(vcf_path: str) -> 'RegionSet':
        regions = []
        with closing(pysam.VariantFile(vcf_path)) as vcf_file:
            for vcf_rec in vcf_file:
                regions.append(Region(chrom=vcf_rec.chrom, start=vcf_rec.start, end=vcf_rec.stop,
                                      orig_start=vcf_rec.start, orig_end=vcf_rec.stop))
                if 'TARGET' in vcf_rec.info and isinstance(vcf_info['TARGET'], int):
                    target_chrom = vcf_rec.info.get('TARGET_CHROM', vcf_rec.chrom)
                    target_start = vcf_rec.info['TARGET'] - 1
                    regions.append(Region(chrom=target_chrom, start=target_start, end=target_start,
                                          orig_start=target_start, orig_end=target_start))
                    
        return RegionSet(regions)

    @staticmethod
    def from_fasta(fasta_path: str, region_kind: str) -> 'RegionSet':
        with pysam.FastaFile(fasta_path) as fasta_file:
            regions = []
            for chrom, chrom_length in zip(fasta_file.references, fasta_file.lengths):
                regions.append(Region(chrom=chrom, start=0, end=chrom_length,
                                      kind=region_kind,
                                      orig_start=0, orig_end=chrom_length))
            return RegionSet(regions)

    def overlaps_region(self, test_region: Region) -> bool:
        chrom_itree = self.chrom2itree[test_region.chrom]
        if (not test_region.is_empty() and
            chrom_itree.overlaps_range(test_region.start, test_region.end)):
            return True
        if (test_region.is_empty() and
            any(iv.begin < test_region.start
                for iv in chrom_itree.at(test_region.start))):
            return True
        chrom_empty_regions = self.chrom2empty_regions[test_region.chrom]
        for reverse in (False, True):
            if any(True for _ in chrom_empty_regions.irange(test_region.head(0),
                                                            test_region.tail(0),
                                                            inclusive=(False, False),
                                                            reverse=reverse)):
                return True
        return False

    def get_region_list(self) -> list[Region]:
        return ([ival.data for chrom_itree in self.chrom2itree.values() for ival in chrom_itree] +
                [empty_region for chrom_empty_regions in self.chrom2empty_regions.values()
                 for empty_region in chrom_empty_regions])

    def filtered(self, region_filter: RegionFilter) -> 'RegionSet':
        """Construct a RegionSet containing regions from self that meet given filter"""

        # TODO: use numpy to filter regions faster
        
        def satisfies_filter(region: Region) -> bool:
            return region_filter.satisfied_for(region)

        return RegionSet(filter(satisfies_filter, self.get_region_list()))

    def add_region_set(self, other: 'RegionSet') -> None:
        for chrom, other_chrom_itree in other.chrom2itree.items():
            self.chrom2itree[chrom].update(other_chrom_itree)
        for chrom, other_chrom_empty_regions in other.chrom2empty_regions.items():
            self.chrom2empty_regions[chrom].update(other_chrom_empty_regions)

    def add_region(self, region: Region) -> None:
        self.add_region_set(RegionSet([region]))

    def chop(self, chop_region: Region) -> list[Region]:
        """Chop off region parts that overlap `chop_region`,
        and return any non-empty remnants created in the process."""

        new_regions: list[Region] = []
        chrom_itree = self.chrom2itree[chop_region.chrom]
        def adjust_region(ival: Interval, is_begin: bool) -> Region:
            assert isinstance(ival.data, Region)
            assert ival.data.chrom == chop_region.chrom
            assert ival.data.start == ival.begin
            assert ival.data.end == ival.end
            if is_begin:
                new_region = ival.data.replace(end=chop_region.start)
            else:
                new_region = ival.data.replace(start=chop_region.end)
            new_regions.append(new_region)
            return new_region
        chrom_itree.chop(chop_region.start, chop_region.end, datafunc=adjust_region)

        if not chop_region.is_empty():
            chrom_empty_regions = self.chrom2empty_regions[chop_region.chrom]
            for region in list(chrom_empty_regions.irange(chop_region.head(0),
                                                          chop_region.tail(0),
                                                          inclusive=(False, False))):
                chrom_empty_regions.discard(region)

        return new_regions

# class class RegionSet

def percent_N(seq: str) -> float:
    return 0 if len(seq) == 0 else (seq.count('N') + seq.count('n')) / len(seq)

def pairwise(iterable: Iterable):
    # pairwise('ABCD') → [('A', 'B'), ('B', 'C'), ('C', 'D')]
    iterator = iter(iterable)
    a = next(iterator, None)
    for b in iterator:
        yield a, b
        a = b

def has_duplicates(it: Iterable) -> bool:
    items = set()
    for item in it:
        if item in items:
            return True
        items.add(item)
    return False

@contextmanager
def error_context(*args):
    context_messages = tuple(map(str, args))
    try:
        yield
    except RuntimeError as e:
        e.args += context_messages
        raise
