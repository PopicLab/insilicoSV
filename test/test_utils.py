from copy import deepcopy
from functools import partial

import pytest

from insilicosv.utils import (
    if_not_none, complement, reverse_complement, divergence,
    generate_seq, pairwise,
    Region, RegionSet, RegionFilter, OverlapMode, n_valid_placements,
    has_duplicates)
                              
def test_if_not_none():
    assert if_not_none(1, 2) == 1
    assert if_not_none(None, 2) == 2
    assert if_not_none(False, 2) == False
    assert if_not_none(True, 2) == True
    assert if_not_none(0, 2) == 0
    assert if_not_none((), 2) == ()
    assert if_not_none([], 2) == []
    assert if_not_none([None], 2) == [None]
    assert if_not_none(None, False) == False
    assert if_not_none(None, True) == True

def test_complement():
    assert complement('') == ''
    assert complement('TCGA') == 'AGCT'

def test_reverse_complement():
    assert reverse_complement('') == ''
    assert reverse_complement('TCGA') == 'TCGA'
    assert reverse_complement('ttt') == 'AAA'
    assert reverse_complement('tNt') == 'ANA'

def test_divergence():
    assert divergence('', 1) == ''
    assert divergence('T', 1) in 'CGA'
    assert divergence('A', 1) in 'CGT'
    assert divergence('G', 1) in 'CTA'
    assert divergence('C', 1) in 'GTA'
    seq = generate_seq(100)
    assert divergence(seq, 0) == seq
    assert divergence('TN', 1) in ('AN', 'GN', 'CN')
    assert divergence('tN', 1) in ('AN', 'GN', 'CN')

def test_pairwise():
    assert list(pairwise('')) == []
    assert list(pairwise('A')) == []
    assert list(pairwise('AB')) == [('A', 'B')]
    assert list(pairwise('ABC')) == [('A', 'B'), ('B', 'C')]

def test_region():
    Reg = partial(Region, chrom="chr21")

    assert Reg(start=0, end=5).length() == 5

    assert Reg(start=5, end=10, kind='ALU') != Reg(start=5, end=10)
    assert Reg(start=5, end=10, kind='ALU', data='3') != Reg(start=5, end=10)

    with pytest.raises(Exception):
        Region(chrom='', start=10, end=15)
    with pytest.raises(Exception):
        Reg(start=10, end=9)

def test_region_set():
    regions1 = [Region('chrA', 1, 5), Region('chrA', 8, 10), Region('chrA', 12, 12),
                Region('chrB', 5, 15)]
    rs1 = RegionSet(deepcopy(regions1))
    assert len(rs1) == len(regions1)
    assert sorted(rs1.get_region_list()) == regions1
    for region in regions1:
        assert region in rs1
        
    assert rs1.overlaps_region(Region('chrA', 1, 5))
    assert rs1.overlaps_region(Region('chrA', 1, 2))
    assert rs1.overlaps_region(Region('chrA', 0, 2))
    assert rs1.overlaps_region(Region('chrA', 0, 8))
    assert rs1.overlaps_region(Region('chrA', 4, 5))
    assert not rs1.overlaps_region(Region('chrC', 4, 5))
    assert not rs1.overlaps_region(Region('chrA', 1, 1))
    assert not rs1.overlaps_region(Region('chrA', 5, 5))
    assert not rs1.overlaps_region(Region('chrA', 12, 12))
    assert not rs1.overlaps_region(Region('chrA', 12, 13))
    assert not rs1.overlaps_region(Region('chrA', 11, 12))
    assert rs1.overlaps_region(Region('chrA', 11, 13))
    assert rs1.overlaps_region(Region('chrA', 2, 2))

    rs1_orig = deepcopy(rs1)
    chop1 = rs1.chop(Region('chrZ', 5, 10))
    assert not chop1
    assert sorted(rs1.get_region_list()) == sorted(rs1_orig.get_region_list())

    chop2 = rs1.chop(Region('chrA', 2, 4))
    assert sorted(chop2) == [Region('chrA', 1, 2), Region('chrA', 4, 5)]
    assert sorted(rs1.get_region_list()) == sorted([
        Region('chrA', 1, 2), Region('chrA', 4, 5), Region('chrA', 8, 10), Region('chrA', 12, 12),
        Region('chrB', 5, 15)])

    chop3 = rs1.chop(Region('chrB', 10, 10))
    assert sorted(chop3) == [Region('chrB', 5, 10), Region('chrB', 10, 15)]
    assert sorted(rs1.get_region_list()) == sorted([
        Region('chrA', 1, 2), Region('chrA', 4, 5), Region('chrA', 8, 10), Region('chrA', 12, 12),
        Region('chrB', 5, 10), Region('chrB', 10, 15)])
    assert not rs1.overlaps_region(Region('chrB', 10, 10))

    chop4 = rs1.chop(Region('chrB', 9, 12))
    assert sorted(chop4) == sorted([Region('chrB', 5, 9), Region('chrB', 12, 15)])
    assert sorted(rs1.get_region_list()) == sorted([
        Region('chrA', 1, 2), Region('chrA', 4, 5), Region('chrA', 8, 10), Region('chrA', 12, 12),
        Region('chrB', 5, 9), Region('chrB', 12, 15)])

    # overlapping regions
    rs1.add_region(Region('chrB', 4, 14))
    assert sorted(rs1.get_region_list()) == sorted([
        Region('chrA', 1, 2), Region('chrA', 4, 5), Region('chrA', 8, 10), Region('chrA', 12, 12),
        Region('chrB', 4, 14), Region('chrB', 5, 9), Region('chrB', 12, 15)])

    rs1.add_region(Region('chrB', 3, 3))
    assert sorted(rs1.get_region_list()) == sorted([
        Region('chrA', 1, 2), Region('chrA', 4, 5), Region('chrA', 8, 10), Region('chrA', 12, 12),
        Region('chrB', 3, 3),
        Region('chrB', 4, 14), Region('chrB', 5, 9), Region('chrB', 12, 15)])

    chop5 = rs1.chop(Region('chrB', 2, 4))
    # we'd expect
    # assert chop5 == []
    # but due to intervaltree issue ( https://github.com/chaimleib/intervaltree/issues/141 )
    # this is currently
    assert chop5 in ([], [Region('chrB', 4, 14)])
    assert sorted(rs1.get_region_list()) == sorted([
        Region('chrA', 1, 2), Region('chrA', 4, 5), Region('chrA', 8, 10), Region('chrA', 12, 12),
        Region('chrB', 4, 14), Region('chrB', 5, 9), Region('chrB', 12, 15)])

    # regions with extra data
    assert Region('chrA', 8, 10) in rs1
    assert Region('chrA', 8, 10, kind='L1') not in rs1
    rs1.add_region_set(RegionSet([Region('chrA', 8, 10, kind='L1')]))
    assert Region('chrA', 8, 10, kind='L1') in rs1
    assert Region('chrA', 8, 10) in rs1
    assert sorted(rs1.get_region_list()) == sorted([
        Region('chrA', 1, 2),
        Region('chrA', 4, 5), Region('chrA', 8, 10), Region('chrA', 8, 10, kind='L1'),
        Region('chrA', 12, 12),
        Region('chrB', 4, 14), Region('chrB', 5, 9), Region('chrB', 12, 15)])

    chop6 = rs1.chop(Region('chrA', 9, 12, kind='ALU'))
    assert sorted(chop6) == sorted([Region('chrA', 8, 9), Region('chrA', 8, 9, kind='L1')])
    assert sorted(rs1.get_region_list()) == sorted([
        Region('chrA', 1, 2),
        Region('chrA', 4, 5), Region('chrA', 8, 9), Region('chrA', 8, 9, kind='L1'),
        Region('chrA', 12, 12),
        Region('chrB', 4, 14), Region('chrB', 5, 9), Region('chrB', 12, 15)])
    

def create_regions():
    return [
        Region(chrom="chr1", start=100, end=200),
        Region(chrom="chr1", start=300, end=400),
        Region(chrom="chr2", start=50, end=150),
        Region(chrom="chr2", start=250, end=350),
        Region(chrom="chr3", start=1000, end=1100)
    ]

def test_region_set_init():
    regions = create_regions()
    region_set = RegionSet(regions)
    
    assert len(region_set) == 5
    assert len(region_set.chrom2itree["chr1"]) == 2
    assert len(region_set.chrom2itree["chr2"]) == 2
    assert len(region_set.chrom2itree["chr3"]) == 1

def test_region_set_contains():
    regions = create_regions()
    region_set = RegionSet(regions)
    
    assert regions[0] in region_set
    assert Region(chrom="chr1", start=100, end=150) not in region_set

def test_region_set_add_region():
    region_set = RegionSet()
    new_region = Region(chrom="chr1", start=100, end=200)
    region_set.add_region(new_region)
    
    assert len(region_set) == 1
    assert new_region in region_set

def test_region_set_add_region_set():
    region_set_1 = RegionSet([Region(chrom="chr1", start=100, end=200)])
    region_set_2 = RegionSet([Region(chrom="chr1", start=300, end=400)])
    
    region_set_1.add_region_set(region_set_2)
    
    assert len(region_set_1) == 2
    assert Region(chrom="chr1", start=300, end=400) in region_set_1

def test_region_set_filtered():
    regions = create_regions()
    region_set = RegionSet(regions)
    
    region_filter = RegionFilter(region_kinds=("L1",), region_length_range=(100, 300))
    filtered_set = region_set.filtered(region_filter)
    
    assert len(filtered_set) == 0  # No regions in `regions` satisfy the filter

    # Add a region that matches the filter
    matching_region = Region(chrom="chr1", start=100, end=300, kind="L1")
    region_set.add_region(matching_region)
    
    filtered_set = region_set.filtered(region_filter)
    assert len(filtered_set) == 1
    assert matching_region in filtered_set

def test_region_set_overlaps_region():
    regions = create_regions()
    region_set = RegionSet(regions)
    
    overlapping_region = Region(chrom="chr1", start=150, end=250)
    non_overlapping_region = Region(chrom="chr1", start=500, end=600)
    
    assert region_set.overlaps_region(overlapping_region) is True
    assert region_set.overlaps_region(non_overlapping_region) is False

def test_region_set_get_region_list():
    regions = create_regions()
    region_set = RegionSet(regions)
    
    region_list = region_set.get_region_list()
    assert len(region_list) == 5
    assert all(isinstance(region, Region) for region in region_list)

def test_region_set_chop():
    regions = [
        Region(chrom="chr1", start=100, end=200),
        Region(chrom="chr1", start=250, end=350),
        Region(chrom="chr1", start=400, end=500)
    ]
    region_set = RegionSet(regions)
    
    chop_region = Region(chrom="chr1", start=150, end=300)
    
    chopped_regions = region_set.chop(chop_region)
    
    # The region from 100 to 200 should be split into two parts (100 to 150) and (300 to 350).
    assert len(chopped_regions) == 2
    assert chopped_regions[0] == Region(chrom="chr1", start=100, end=150)
    assert chopped_regions[1] == Region(chrom="chr1", start=300, end=350)
    
    # The original region set should reflect the changes.
    assert len(region_set) == 3  # 100-150, 300-350, and 400-500 remain

def test_region_set_add_empty_region():
    empty_region = Region(chrom="chr1", start=100, end=100)  # Empty region
    region_set = RegionSet([empty_region])
    
    assert empty_region in region_set
    assert len(region_set.chrom2empty_regions["chr1"]) == 1
    assert len(region_set) == 1

def test_region_set_assert_valid():
    valid_regions = create_regions()
    region_set = RegionSet(valid_regions)
    
    region_set.assert_valid()

def test_region_set_empty():
    region_set = RegionSet([])
    assert len(region_set) == 0
    assert region_set.get_region_list() == []

def test_region_set_large_region_handling():
    large_region = Region(chrom="chr1", start=1000, end=100000)
    region_set = RegionSet([large_region])
    
    assert len(region_set) == 1
    assert large_region in region_set

def test_region_set_from_beds(tmpdir):
    bed_file = tmpdir.join("test.bed")
    bed_file.write("chr1\t100\t200\tkind1\nchr1\t300\t400\tkind2\n")

    region_set = RegionSet.from_beds([str(bed_file)])

    assert len(region_set) == 2
    assert Region(chrom="chr1", start=100, end=200, kind="kind1", orig_start=100, orig_end=200) in region_set
    assert Region(chrom="chr1", start=300, end=400, kind="kind2", orig_start=300, orig_end=400) in region_set

def test_region_set_from_vcf(tmpdir):
    vcf_file = tmpdir.join("test.vcf")
    vcf_file.write("##fileformat=VCFv4.2\n"
                   "##contig=<ID=chr1,length=1000>\n"
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "chr1\t100\t.\tA\tT\t.\t.\t.\n"
                   "chr1\t300\t.\tG\tC\t.\t.\t.\n")
    
    region_set = RegionSet.from_vcf(str(vcf_file))
    
    assert len(region_set) == 2, f'{region_set.get_region_list()=}'
    assert Region(chrom="chr1", start=99, end=100, orig_start=99, orig_end=100) in region_set
    assert Region(chrom="chr1", start=299, end=300, orig_start=299, orig_end=300) in region_set

def test_region_set_from_fasta(tmpdir):
    fasta_file = tmpdir.join("test.fasta")
    fasta_file.write(">chr1\n" + "A" * 100 + "\n>chr2\n" + "C" * 200 + "\n")
    
    region_set = RegionSet.from_fasta(str(fasta_file), "genome")
    
    assert len(region_set) == 2

    assert Region(chrom="chr1", start=0, end=100, kind="genome", orig_start=0, orig_end=100) in region_set
    assert Region(chrom="chr2", start=0, end=200, kind="genome", orig_start=0, orig_end=200) in region_set

def test_n_valid_placements():
    assert n_valid_placements(overlap_mode=OverlapMode.CONTAINED, region_length=1, anchor_length=1) == 1
    assert n_valid_placements(overlap_mode=OverlapMode.CONTAINED, region_length=1, anchor_length=2) == 0
    assert n_valid_placements(overlap_mode=OverlapMode.CONTAINED, region_length=3, anchor_length=2) == 2
    assert n_valid_placements(overlap_mode=OverlapMode.CONTAINED, region_length=1, anchor_length=0) == 0
    assert n_valid_placements(overlap_mode=OverlapMode.CONTAINED, region_length=2, anchor_length=0) == 1
    assert n_valid_placements(overlap_mode=OverlapMode.CONTAINED, region_length=3, anchor_length=0) == 2

    assert n_valid_placements(overlap_mode=OverlapMode.PARTIAL, region_length=1, anchor_length=2) == 0
    assert n_valid_placements(overlap_mode=OverlapMode.PARTIAL, region_length=2, anchor_length=2) == 2
    assert n_valid_placements(overlap_mode=OverlapMode.PARTIAL, region_length=3, anchor_length=1) == 0
    assert n_valid_placements(overlap_mode=OverlapMode.PARTIAL, region_length=3, anchor_length=2) == 2
    assert n_valid_placements(overlap_mode=OverlapMode.PARTIAL, region_length=3, anchor_length=3) == 4

def test_has_duplicates():
    assert not has_duplicates([])
    assert not has_duplicates([1])
    assert not has_duplicates([1,2])
    assert has_duplicates([1,2,1])
    assert has_duplicates('aa')
    assert not has_duplicates(())
    assert has_duplicates([(), 1, 'a', ()])
