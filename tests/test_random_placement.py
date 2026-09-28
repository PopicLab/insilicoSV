import random
from types import SimpleNamespace
import yaml

import pytest

from insilicosv.simulate import SVSimulator
from insilicosv.utils import OverlapMode, Region, RegionSet

DUMMY_VSET_CONFIG = {"config_descr": "test_sv_config", "intersv_dist_min": 0, "max_tries": 50, "random_seed": 0, "allow_hap_overlap": False}

def _make_regions(small, large, num_regions):
    small = [small for _ in range(num_regions // 2)]
    large = [large for _ in range(num_regions // 2)]
    return small + large

def _run_overlap_mode_test(tmp_path, overlap_mode, num_svs, num_regions, min_ratio, max_ratio):
    small_region_length = 1000
    large_region_length = 10000
    regions = _make_regions(small_region_length, large_region_length, num_regions)
    dist_regions = 100
    anchor_length = 50

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
    for region_len in regions:
        rois.append(Region(
            "chr1", 
            start_x, 
            start_x + region_len, 
            orig_start=start_x, 
            orig_end=start_x + region_len, 
            region_type=str(region_len)
        ))
        start_x += region_len + dist_regions

    if overlap_mode == OverlapMode.CONTAINED:
        simulator.rois_overlap = {0: RegionSet(rois)}
        simulator.union_rois_overlap = {0: RegionSet()}
        simulator.union_rois_overlap[0].build_union_tree(simulator.rois_overlap[0])
    else:
        simulator.rois_overlap = {0: rois}
        random.seed(0)
        random.shuffle(simulator.rois_overlap[0])

    roi_filter = SimpleNamespace(
        region_length_range=(anchor_length, None),
        satisfied_for=lambda region: True,
    )

    random.seed(0)
    
    small_count = 0
    large_count = 0
    roi_index = 0

    for _ in range(num_svs):
        overlap_roi, ref_roi, roi_index, random_position = simulator.get_overlap_region(
            sv_category=0,
            roi_index=roi_index,
            reference_regions=simulator.reference_regions,
            hap_id=0,
            anchor_length=anchor_length if overlap_mode == OverlapMode.CONTAINED else None,
            overlap_mode=overlap_mode,
            roi_filter=roi_filter,
        )

        assert overlap_roi is not None
        assert ref_roi is not None

        anchor_len_to_pass = anchor_length if overlap_mode == OverlapMode.CONTAINED else overlap_roi.length()

        anchor_region, anchor_ref = simulator.choose_anchor_placement(
            roi=overlap_roi,
            ref_roi=ref_roi,
            anchor_length=anchor_len_to_pass,
            overlap_mode=overlap_mode,
            region_length_range=roi_filter.region_length_range,
            blacklist_regions=None,
            random_position=random_position,
        )

        assert anchor_region is not None
        assert anchor_ref is ref_roi
        assert anchor_region.chrom == overlap_roi.chrom
        
        if overlap_mode == OverlapMode.CONTAINED:
            assert anchor_region.length() == anchor_length
            assert anchor_region.start >= overlap_roi.start
            assert anchor_region.end <= overlap_roi.end
        elif overlap_mode == OverlapMode.EXACT:
            assert anchor_region.start == overlap_roi.start
            assert anchor_region.end == overlap_roi.end
        
        # Tally sizes based on the region's tagged original length
        if overlap_roi.region_type == str(small_region_length):
            small_count += 1
        elif overlap_roi.region_type == str(large_region_length):
            large_count += 1
        else:
            pytest.fail(f"Unexpected region size selected: {overlap_roi.region_type}")

    assert small_count + large_count == num_svs
    
    # Assert specific distributions based on how the overlap mode selects regions
    assert small_count > 0, "No SVs were placed in small regions"
    ratio = large_count / small_count
    assert min_ratio <= ratio <= max_ratio, f"Expected ratio between {min_ratio} and {max_ratio}, got {ratio:.2f} ({large_count} large / {small_count} small)"
    print(f"Overlap mode {overlap_mode.name}: {small_count} small, {large_count} large, ratio = {ratio:.2f}")


def test_random_placement_contained(tmp_path):
    # CONTAINED samples uniformly by length across the union tree, so 1000bp vs 100bp regions = ~10:1 ratio. 
    # Run 1000 times to get a stable statistical distribution.
    _run_overlap_mode_test(
        tmp_path, 
        overlap_mode=OverlapMode.CONTAINED, 
        num_svs=10000, 
        num_regions=100,
        min_ratio=8.0, 
        max_ratio=12.0
    )


def test_random_placement_exact(tmp_path):
    # EXACT iterates through a shuffled flat list of the 100 ROIs directly, so selection is strictly 1:1.
    # Run 100 times to exactly consume the 100 ROIs (50 small, 50 large).
    _run_overlap_mode_test(
        tmp_path, 
        overlap_mode=OverlapMode.EXACT, 
        num_svs=1000,
        num_regions=10000, 
        min_ratio=0.7, 
        max_ratio=1.3
    )