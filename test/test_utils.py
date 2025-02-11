from copy import deepcopy
from functools import partial

import pytest

from insilicosv.utils import (
    if_not_none, complement, reverse_complement, divergence,
    generate_seq, pairwise,
    Region, RegionSet, RegionFilter, OverlapMode,
    has_duplicates)

def region_length(region):
    return sum(len(chrom_itree) for chrom_itree in region.chrom2itree.values())

def test_chop():
    regions1 = [Region('chrA', 1, 5), Region('chrA', 8, 10),
                Region('chrB', 5, 15)]
    rs1 = RegionSet(deepcopy(regions1))
    rs1.chop(Region('chrA', 1, 5))
    assert len(rs1.chrom2itree['chrA']) == 1
    assert len(rs1.chrom2itree['chrB']) == 1
    rs1.chop(Region('chrB', 5, 15))
    assert len(rs1.chrom2itree['chrB']) == 0
    print('final tree', rs1.chrom2itree)


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
    regions1 = [Region('chrA', 1, 5), Region('chrA', 8, 10), Region('chrA', 12, 13),
                Region('chrB', 5, 15)]
    rs1 = RegionSet(deepcopy(regions1))
    assert region_length(rs1) == len(regions1)
    assert sorted(rs1.get_region_list()) == regions1
    for region in regions1:
        assert region in rs1

    rs1_orig = deepcopy(rs1)
    rs1.chop(Region('chrZ', 5, 10))

    rs1.chop(Region('chrA', 2, 4))


    rs1.chop(Region('chrB', 10, 10))


    rs1.chop(Region('chrB', 9, 12))


    # overlapping regions
    rs1.add_region(Region('chrB', 4, 14))


    rs1.chop(Region('chrB', 2, 4))
    # we'd expect
    # assert chop5 == []
    # but due to intervaltree issue ( https://github.com/chaimleib/intervaltree/issues/141 )
    # this is currently

    # regions with extra data
    assert Region('chrA', 8, 10) in rs1
    rs1.add_region_set(RegionSet([Region('chrA', 8, 10, kind='L1')]))
    assert Region('chrA', 8, 10, kind='L1') in rs1
    assert Region('chrA', 8, 10) in rs1

    rs1.chop(Region('chrA', 9, 12, kind='ALU'))

    

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
    
    assert region_length(region_set) == 5
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
    
    assert region_length(region_set)== 1
    assert new_region in region_set

def test_region_set_add_region_set():
    region_set_1 = RegionSet([Region(chrom="chr1", start=100, end=200)])
    region_set_2 = RegionSet([Region(chrom="chr1", start=300, end=400)])
    
    region_set_1.add_region_set(region_set_2)
    
    assert region_length(region_set_1) == 2
    assert Region(chrom="chr1", start=300, end=400) in region_set_1

def test_region_set_filtered():
    regions = create_regions()
    region_set = RegionSet(regions)
    
    region_filter = RegionFilter(region_kinds=("L1",), region_length_range=(100, 300))
    filtered_set = region_set.filtered(region_filter)
    
    assert region_length(filtered_set)== 0  # No regions in `regions` satisfy the filter

    # Add a region that matches the filter
    matching_region = Region(chrom="chr1", start=100, end=300, kind="L1")
    region_set.add_region(matching_region)
    
    filtered_set = region_set.filtered(region_filter)
    assert region_length(filtered_set)== 1
    assert matching_region in filtered_set

def test_region_set_get_region_list():
    regions = create_regions()
    region_set = RegionSet(regions)
    
    region_list = region_set.get_region_list()
    assert len(region_list) == 5
    assert all(isinstance(region, Region) for region in region_list)

def test_region_set_empty():
    region_set = RegionSet([])
    assert region_length(region_set) == 0
    assert region_set.get_region_list() == []

def test_region_set_large_region_handling():
    large_region = Region(chrom="chr1", start=1000, end=100000)
    region_set = RegionSet([large_region])
    
    assert region_length(region_set) == 1
    assert large_region in region_set

def test_region_set_from_beds(tmpdir):
    bed_file = tmpdir.join("test.bed")
    bed_file.write("chr1\t100\t200\tkind1\nchr1\t300\t400\tkind2\n")

    region_set = RegionSet.from_beds([str(bed_file)], True)

    assert region_length(region_set) == 2
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
    
    assert region_length(region_set) == 2, f'{region_set.get_region_list()=}'
    assert Region(chrom="chr1", start=99, end=100, orig_start=99, orig_end=100) in region_set
    assert Region(chrom="chr1", start=299, end=300, orig_start=299, orig_end=300) in region_set

def test_region_set_from_fasta(tmpdir):
    fasta_file = tmpdir.join("test.fasta")
    fasta_file.write(">chr1\n" + "A" * 100 + "\n>chr2\n" + "C" * 200 + "\n")
    
    region_set = RegionSet.from_fasta(str(fasta_file), "genome")
    
    assert region_length(region_set) == 2

    assert Region(chrom="chr1", start=0, end=100, kind="genome", orig_start=0, orig_end=100) in region_set
    assert Region(chrom="chr2", start=0, end=200, kind="genome", orig_start=0, orig_end=200) in region_set

def test_has_duplicates():
    assert not has_duplicates([])
    assert not has_duplicates([1])
    assert not has_duplicates([1,2])
    assert has_duplicates([1,2,1])
    assert has_duplicates('aa')
    assert not has_duplicates(())
    assert has_duplicates([(), 1, 'a', ()])
