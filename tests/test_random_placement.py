import random
from types import SimpleNamespace
import yaml

import pytest

from insilicosv.simulate import SVSimulator
from insilicosv.utils import OverlapMode, Region, RegionSet

DUMMY_VSET_CONFIG = {"config_descr": "test_sv_config", "intersv_dist_min": 0, "max_tries": 50, "random_seed": 0, "allow_hap_overlap": False}

def _make_regions(small, large):
    small = [small for _ in range(50)]
    large = [large for _ in range(50)]
    return small + large

def test_simulator_get_overlap_region_and_choose_anchor_placement(tmp_path):
    small_region_length = 1000
    large_region_length = 10000
    regions = _make_regions(small_region_length, large_region_length)
    dist_regions = 100

    fasta_path = tmp_path / "ref.fa"
    fasta_path.write_text(">chr1\n" + "A" * sum(regions) + "A" * (dist_regions * len(regions)) + "\n")

    config_path = tmp_path / "config.yaml"
    
    config_data = {
        "reference": str(fasta_path),
        "variant_sets": [DUMMY_VSET_CONFIG],
        "max_tries": 50,
        "random_seed": 0
    }
    config_path.write_text(yaml.dump(config_data))

    simulator = SVSimulator(str(config_path))
    
    simulator.reference_regions = RegionSet.from_fasta(
        str(fasta_path), 
        0, 
        region_region_type="_reference_",
        allow_hap_overlap=simulator.allow_hap_overlap
    )

    rois = []
    start_x = 0
    for region in regions:
        rois.append(Region("chr1", start_x, start_x + region, region_type=str(region)))
        start_x += region + dist_regions

    simulator.rois_overlap = {0: RegionSet(rois)}
    simulator.union_rois_overlap = {0: RegionSet()}
    simulator.union_rois_overlap[0].build_union_tree(simulator.rois_overlap[0])

    roi_filter = SimpleNamespace(
        region_length_range=(50, None),
        satisfied_for=lambda region: True,
    )

    random.seed(0)
    
    num_svs = 10000
    small_count = 0
    large_count = 0

    for _ in range(num_svs):
        overlap_roi, ref_roi, _, random_position = simulator.get_overlap_region(
            sv_category=0,
            roi_index=0,
            reference_regions=simulator.reference_regions,
            hap_id=0,
            anchor_length=50,
            overlap_mode=OverlapMode.CONTAINED,
            roi_filter=roi_filter,
        )

        assert overlap_roi is not None
        assert ref_roi is not None

        anchor_region, anchor_ref = simulator.choose_anchor_placement(
            roi=overlap_roi,
            ref_roi=ref_roi,
            anchor_length=50,
            overlap_mode=OverlapMode.CONTAINED,
            region_length_range=roi_filter.region_length_range,
            blacklist_regions=None,
            random_position=random_position,
        )

        assert anchor_region is not None
        assert anchor_ref is ref_roi
        assert anchor_region.length() == 50
        assert anchor_region.start >= overlap_roi.start
        assert anchor_region.end <= overlap_roi.end
        assert anchor_region.chrom == overlap_roi.chrom
        
        if overlap_roi.region_type == str(small_region_length):
            small_count += 1
        elif overlap_roi.region_type == str(large_region_length):
            large_count += 1
        else:
            pytest.fail(f"Unexpected region size selected: {overlap_roi.region_type}")

    assert small_count + large_count == num_svs
    
    assert small_count > 0, "No SVs were placed in small regions"
    ratio = large_count / small_count
    assert 8.0 <= ratio <= 12.0, f"Expected roughly 10x SVs in large regions, got ratio {ratio:.2f} ({large_count} large / {small_count} small)"
    print(f"SV placement counts: {small_count} small, {large_count} large, ratio {ratio:.2f}")