import os
import shutil
import sys
import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Any, Union

import numpy as np
import pytest
import yaml
from pysam import FastaFile

from insilicosv import utils
from insilicosv.simulate import SVSimulator
from insilicosv.sv_defs import SV, Transform, TransformType
from insilicosv.utils import Region
import logging
logging.getLogger("insilicosv.output").setLevel(logging.ERROR)

# ==========================================
# FILE PATH CONSTANTS
# ==========================================
OVERLAP_BED_1 = "tests/inputs/example_overlap_events.bed"
OVERLAP_BED_2 = "tests/inputs/example_overlap_events_2.bed"
OVERLAP_BED_5 = "tests/inputs/example_overlap_events_5.bed"
OVERLAP_BED_7 = "tests/inputs/example_overlap_events_7.bed"
OVERLAP_BED_8 = "tests/inputs/example_overlap_events_8.bed"
OVERLAP_BED_9 = "tests/inputs/example_overlap_events_9.bed"
OVERLAP_BED_10 = "tests/inputs/example_overlap_events_10.bed"
OVERLAP_BED_11 = "tests/inputs/example_overlap_events_11.bed"
OVERLAP_BED_12 = "tests/inputs/example_overlap_events_12.bed"
OVERLAP_BED_13 = "tests/inputs/example_overlap_events_13.bed"
OVERLAP_BED_14 = "tests/inputs/example_overlap_events_14.bed"
OVERLAP_BED_15 = "tests/inputs/example_overlap_events_15.bed"
OVERLAP_VCF = "tests/inputs/example_overlap.vcf"
MOTIF_BED = "tests/inputs/motif.bed"

IMPORT_DEL = "tests/inputs/import_del.vcf"
IMPORT_SNP = "tests/inputs/import_snp.vcf"
IMPORT_INV = "tests/inputs/import_inv.vcf"
IMPORT_TEST = "tests/inputs/import_test.vcf"
IMPORT_SNP_OVERLAP = "tests/inputs/import_snp_overlap.vcf"

EXCLUDE_BED = "tests/inputs/exclude.bed"
INSERTION_ORDER_EXACT = "tests/inputs/test_insertion_order_exact.bed"
AVOID_INTERVAL = "tests/inputs/example_avoid_interval.vcf"
AVOID_INTERVAL_2 = "tests/inputs/example_avoid_interval_2.bed"
AVOID_INTERVAL_3 = "tests/inputs/example_avoid_interval_3.bed"


# ==========================================
# DATA CLASSES & HELPERS
# ==========================================

@dataclass
class SVTestCase:
    """Pure data object representing a single configuration to test."""
    id: str
    ref_seqs: Dict[str, str]
    variant_sets: List[Dict]
    expected_outputs: List[str] = field(default_factory=list)
    config_overrides: Dict[str, Any] = field(default_factory=dict)
    heterozygous: bool = False


def make_case(cid: str, ref: Dict[str, str], vsets: List[Dict], expected: List[str] = None, hetero=False, **kwargs):
    """Helper to cleanly build standard SVTestCases."""
    return SVTestCase(cid, ref, vsets, expected or [], kwargs, hetero)


def make_rt_case(cid: str, ref: Union[str, List[str]], vsets: List[Dict], expected: List[str] = None, hetero=False, **kwargs):
    """Helper that exactly mimics the config generation of the original run_test() method."""
    if isinstance(ref, str):
        ref_dict = {"chrTest0": ref}
    else:
        ref_dict = {f"chrTest{i}": seq for i, seq in enumerate(ref)}
    
    overrides = {"homozygous_only": True, "min_intersv_dist": 0, "random_seed": 55}
    overrides.update(kwargs)
    
    # Automatically add number=1 to variant sets if missing
    for vset in vsets:
        if "number" not in vset and "import" not in vset:
            vset["number"] = 1
            
    return SVTestCase(cid, ref_dict, vsets, expected or [], overrides, hetero)


def is_overlapping(event_ranges, addition, strictly_partial=False):
    for event in event_ranges:
        if event[1] > addition[0] and event[0] < addition[1]:
            return True if not strictly_partial else event != addition
    return False


def get_span(sv: SV) -> Region:
    assert sv.is_placed()
    regions = sorted(sv.get_regions())
    chrom = regions[0].chrom
    assert all(region.chrom == chrom for region in regions)
    return Region(chrom=chrom, start=regions[0].start, end=regions[-1].end)


def sv_source_segments(sv):
    return sorted([(op.source_region.start, op.source_region.end) for op in sv.operations])


# ==========================================
# EXECUTION ENGINES
# ==========================================

def run_seed_search(case: SVTestCase, tmp_path, allow_fail=False):
    """
    Tries different random seeds until ALL expected outputs are generated at least once.
    Fails safely if unexpected outputs are generated or if it loops infinitely.
    """
    results_seen = set()
    expected = set(case.expected_outputs)
    attempt_num = 0
    max_attempts = max(len(expected) * 200, 200) # Enough iterations to find all branches
    current_seed = case.config_overrides.get("random_seed", 55)

    while not expected.issubset(results_seen) and attempt_num < max_attempts:
        current_seed += 15
        attempt_num += 1

        run_dir = tmp_path / f"run_{attempt_num}"
        run_dir.mkdir()
        
        ref_file = run_dir / "ref.fna"
        par_file = run_dir / "par.yaml"
        hap1 = run_dir / "sim.hapA.fa"
        hap2 = run_dir / "sim.hapB.fa"

        with open(ref_file, "w") as f:
            for chrom, seq in case.ref_seqs.items():
                f.write(f">{chrom}\n{seq}\n")
                
        with FastaFile(str(ref_file)): # Opens, creates the index, and cleanly closes
            pass

        config = {
            "reference": str(ref_file),
            "variant_sets": case.variant_sets
        }
        config.update(case.config_overrides)
        config["random_seed"] = current_seed  

        #Convert to absolute paths
        for key in ["overlap_regions", "blacklist_regions", "novel_insertions"]:
            if key in config:
                if isinstance(config[key], list):
                    config[key] = [os.path.abspath(p) for p in config[key]]
                elif isinstance(config[key], str):
                    config[key] = os.path.abspath(config[key])
                    
        for vset in config.get("variant_sets", []):
            for key in ["import", "novel_insertions"]:
                if key in vset:
                    if isinstance(vset[key], list):
                        vset[key] = [os.path.abspath(p) for p in vset[key]]
                    elif isinstance(vset[key], str):
                        vset[key] = os.path.abspath(vset[key])

        with open(par_file, "w") as f:
            yaml.dump(config, f, default_flow_style=False)

        try:
            sim = SVSimulator(config_path=str(par_file))
            sim.run()

            with FastaFile(str(hap1)) as f1, FastaFile(str(hap2)) as f2:
                frags1 = [f1.fetch(r) for r in f1.references]
                frags2 = [f2.fetch(r) for r in f2.references]

            frag1 = frags1[0] if len(frags1) == 1 else tuple(frags1)
            frag2 = frags2[0] if len(frags2) == 1 else tuple(frags2)

            if case.heterozygous:
                if frag1 in expected: results_seen.add(frag1)
                if frag2 in expected: results_seen.add(frag2)
                if frag1 not in expected and frag2 not in expected:
                    results_seen.update([frag1, frag2])
            else:
                results_seen.update([frag1, frag2])

        except Exception as e:
            if not allow_fail:
                raise e
            
        # Optimization: If no expected outputs are defined, just run once for crash-testing
        if not expected:
            return

    missing = expected - results_seen
    assert not missing, f"[{case.id}] Missing outputs after {attempt_num} tries: {missing}"
    
    unexpected = results_seen - expected
    assert not unexpected, f"[{case.id}] Found unexpected outputs: {unexpected}"


def run_single_sv_case(case: SVTestCase, tmp_path):
    """Runs a configuration exactly once (for internal inspection tests)."""
    ref_file = tmp_path / "ref.fna"
    par_file = tmp_path / "par.yaml"
    hap1 = tmp_path / "sim.hapA.fa"
    hap2 = tmp_path / "sim.hapB.fa"

    with open(ref_file, "w") as f:
            for chrom, seq in case.ref_seqs.items():
                f.write(f">{chrom}\n{seq}\n")
                
    with FastaFile(str(ref_file)): # Opens, creates the index, and cleanly closes
        pass

    config = {
        "reference": str(ref_file),
        "random_seed": case.config_overrides.get("random_seed", 88),
        "variant_sets": case.variant_sets
    }
    config.update(case.config_overrides)

    #Convert to absolute paths
    for key in ["overlap_regions", "blacklist_regions", "novel_insertions"]:
        if key in config:
            if isinstance(config[key], list):
                config[key] = [os.path.abspath(p) for p in config[key]]
            elif isinstance(config[key], str):
                config[key] = os.path.abspath(config[key])
                
    for vset in config.get("variant_sets", []):
        for key in ["import", "novel_insertions"]:
            if key in vset:
                if isinstance(vset[key], list):
                    vset[key] = [os.path.abspath(p) for p in vset[key]]
                elif isinstance(vset[key], str):
                    vset[key] = os.path.abspath(vset[key])

    with open(par_file, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    sim = SVSimulator(config_path=str(par_file))
    sim.run()

    with FastaFile(str(hap1)) as f1, FastaFile(str(hap2)) as f2:
        frags1 = [f1.fetch(r) for r in f1.references]
        frags2 = [f2.fetch(r) for r in f2.references]

    frag1 = frags1[0] if len(frags1) == 1 else tuple(frags1)
    frag2 = frags2[0] if len(frags2) == 1 else tuple(frags2)

    return frag1, frag2, sim.svs


# ==========================================
# TEST DATA DEFINITIONS
# ==========================================

SNPS_DATA = [
    make_case("snp_1", {"chr21": "CA"}, [{"type": "SNP", "number": 1}],
              ['AA', 'GA', 'TA', 'CC', 'CG', 'CT'], homozygous_only=True),
    make_case("snp_multi", {"chr21": "CTGTTGACCG"}, [{"type": "SNP", "number": 4}],
              homozygous_only=True) 
]

SIMPLE_DELS_DATA = [
    make_case("del_13", {"Chromosome19": "CACTATCTCTCCGAT"}, [{"type": "DEL", "number": 1, "length_ranges": [[13, 13]]}], ['CA', 'CT', 'AT'], homozygous_only=True),
    make_case("del_14", {"Chromosome19": "CACTATCTCTCCGAT"}, [{"type": "DEL", "number": 1, "length_ranges": [[14, 14]]}], ['C', 'T'], homozygous_only=True)
]

SIMPLE_DUPS_DATA = [
    make_case("dup_2", {"Chromosome19": "CA"}, [{"type": "DUP", "number": 1, "length_ranges": [[2, 2]]}], ['CACA'], homozygous_only=True),
    make_case("dup_2_cat", {"Chromosome19": "CAT"}, [{"type": "DUP", "number": 1, "length_ranges": [[2, 2]]}], ['CACAT', 'CATAT'], homozygous_only=True),
    make_case("dup_1", {"Chromosome19": "C"}, [{"type": "DUP", "number": 1, "length_ranges": [[1, 1]]}], ['CC'], homozygous_only=True)
]

SIMPLE_INSS_DATA = [
    make_case("ins_1", {"Chromosome19": "CA"}, [{"type": "INS", "number": 1, "length_ranges": [[5, 5]]}], homozygous_only=True)
]

SIMPLE_INVS_DATA = [
    make_case("inv_2", {"Chromosome19": "CA"}, [{"type": "INV", "number": 1, "length_ranges": [[2, 2]]}], ['TG'], homozygous_only=True),
    make_case("inv_1", {"Chromosome19": "C"}, [{"type": "INV", "number": 1, "length_ranges": [[1, 1]]}], ['G'], homozygous_only=True)
]

MULTI_INS_DATA = [
    make_case("multi_ins_1", {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}, 
              [{"type": "INS", "number": 1, "length_ranges": [[5, 5]]}, 
               {"type": "delINV", "number": 1, "length_ranges": [[5, 5]]}, 
               {"type": "INS", "number": 1, "length_ranges": [[5, 5]]}], homozygous_only=True)
]

OVERLAP_SIMPLE_DATA = [
    make_case("os_0", {"chr21": "CTCCGTCGTA"}, [{"type": "DEL", "number": 1, "length_ranges": [[None, None]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact"}], overlap_regions=[OVERLAP_BED_1]),
    make_case("os_1", {"chr21": "CTCCGTCGTA"}, [{"type": "DUP", "number": 1, "length_ranges": [[None, None]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact"}], overlap_regions=[OVERLAP_BED_1]),
    make_case("os_2", {"chr21": "CTCCGTCGTACTAAGTCGTA"}, [{"type": "DEL", "number": 1, "length_ranges": [[None, None]], "overlap_region_length_range": [2, 5], "overlap_mode": "exact", "overlap_region_type": [["L1PA15"], ["L1PA15"]]}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2])
]

COMPLEX_OVERLAP_DATA = [
    make_case("cplx_0", {"chr21": "CTGAT"}, [{"type": "(A)_->A_A", "number": 1, "length_ranges": [[None, None], [1, 1]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact", "overlap_region_type": [["L1HS"], ['L1HS']]}], ['CTGATGA', 'CGATGAT'], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2], hetero=True),
    make_case("cplx_1", {"chr21": "CTGATATGGAC"}, [{"type": "(A)_->_A", "number": 1, "length_ranges": [[None, None], [1, 1]], "overlap_mode": "exact", "overlap_region_length_range": [4, 6], "overlap_region_type": [["L1HS"], ['L1HS']]}, {"type": "(A)_->A_a", "number": 1, "length_ranges": [[None, None], [1, 1]], "overlap_mode": "exact", "overlap_region_length_range": [1, 1], "overlap_region_type": [["AluSz6"], ['AluSz6']]}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2], hetero=True),
    make_case("cplx_2", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "(A)BC->b", "number": 1, "length_ranges": [[None, None], [3, 3], [3, 3]], "overlap_mode": "exact", "overlap_region_length_range": [3, 3]}], overlap_regions=[OVERLAP_BED_2], hetero=True),
    make_case("cplx_3", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "(A)B->b", "number": 1, "length_ranges": [[None, None], [3, 3]], "overlap_mode": "exact", "overlap_region_length_range": [3, 3], "overlap_region_type": "ALR"}, {"type": "(A)B->a", "number": 1, "length_ranges": [[None, None], [2, 2]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact"}], overlap_regions=[OVERLAP_BED_10], hetero=True),
    make_case("cplx_4", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "(A)BC->b", "number": 1, "length_ranges": [[3, 3], [3, 3], [3, 3]], "overlap_mode": "partial"}], overlap_regions=[OVERLAP_BED_13], hetero=True),
    make_case("cplx_5", {"chr21": "CCTGATCTGATCTGATCTGATCTGATTGAT"}, [{"type": "A_()->A_A", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_region_type": [["L1PA15"], ["L1PA15"]], "overlap_mode": "contained"}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2], hetero=True),
    make_case("cplx_6", {"chr21": "CCTGATCTGATCTGATCTGATCTGATTGAT"}, [{"type": "A_()->A_a", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_region_type": [["L1PA15"], ["L1PA15"]], "overlap_mode": "contained"}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2], hetero=True),
    make_case("cplx_7", {"chr21": "CCTGATCTGATCTGATCTGATCTGATTGAT"}, [{"type": "A(_)->_A", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_region_type": [["L1PA15"], ["L1PA15"]], "overlap_mode": "contained"}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2], hetero=True),
    make_case("cplx_8", {"chr21": "CCTGATATGGACCTGATATGGACTGATATGGAC"}, [{"type": "A_()->A_a", "number": 1, "length_ranges": [[2, 2], [3, 3]], "overlap_region_type": [["ALR"], ["ALR"]], "overlap_mode": "contained"}, {"type": "nrTRA", "number": 1, "length_ranges": [[4, 6], [1, 1]]}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2], max_tries=1000, hetero=True),
    make_case("cplx_9", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "delINVdel", "number": 1, "length_ranges": [[2, 2], [3, 3], [1, 1]], "overlap_mode": "contained", "overlap_region_type": ["L1PA15"]}, {"type": "DEL", "number": 1, "length_ranges": [[2, 6]]}], overlap_regions=[OVERLAP_BED_11], max_tries=1000, hetero=True),
    make_case("cplx_10", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "delINV", "number": 1, "length_ranges": [[3, 3], [3, 3]], "overlap_mode": "contained", "overlap_region_type": ["L1PA15"]}], overlap_regions=[OVERLAP_BED_11], max_tries=1000, hetero=True),
    make_case("cplx_11", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "delINVdel", "number": 1, "length_ranges": [[2, 6], [3, 5], [2, 9]], "overlap_mode": "partial", "overlap_region_type": ["L1PA15"]}], overlap_regions=[OVERLAP_BED_11], hetero=True),
    make_case("cplx_12", {"chr21": "GGACCTATTAG"}, [{"type": "A(BC)__ -> A_C_B", "number": 1, "length_ranges": [[1, 1], [2, 2], [1, 1], [1, 1], [1, 1]], "overlap_mode": "partial", "overlap_region_type": ["L1P"]}], ["GGTCAACTTAG"], overlap_regions=[OVERLAP_BED_14], homozygous_only=True),
    make_case("cplx_13", {"chr21": "GGACCTATTAG"}, [{"type": "()_ABC_ -> A_C_B", "number": 1, "length_ranges": [[1, 1], [1, 1], [2, 2], [1, 1], [1, 1]], "overlap_mode": "contained", "overlap_region_type": ["L1P"]}], ["GAGTACCTTAG", "GGCAATCTTAG"], overlap_regions=[OVERLAP_BED_14], homozygous_only=True),
    make_case("cplx_14", {"chr21": "GGACCTATTAG"}, [{"type": "(_)ABC_ -> A_C_B", "number": 1, "length_ranges": [[2, 2], [1, 1], [2, 2], [1, 1], [1, 1]], "overlap_mode": "partial", "overlap_region_type": ["L1P"]}], ["AGGTACCTTAG", "GGCACTTTAAG"], overlap_regions=[OVERLAP_BED_14], homozygous_only=True),
]

DISPERSION_DATA = [
    make_case("nrTRA", {"Chromosome19": "CT"}, [{"type": "nrTRA", "number": 1, "length_ranges": [[1, 1], [1, 1]]}], ['TC'], homozygous_only=True),
    make_case("dDUP", {"Chromosome19": "CT"}, [{"type": "dDUP", "number": 1, "length_ranges": [[1, 1], [1, 1]]}], ['CTC', 'TCT'], homozygous_only=True),
    make_case("INV_dDUP", {"Chromosome19": "CT"}, [{"type": "INV_dDUP", "number": 1, "length_ranges": [[1, 1], [1, 1]]}], ['CTG', 'ACT'], homozygous_only=True),
    make_case("dDUP_iDEL", {"Chromosome19": "CTG"}, [{"type": "dDUP_iDEL", "number": 1, "length_ranges": [[1, 1], [1, 1], [1, 1]]}], ['CTC', 'GTG'], homozygous_only=True),
    make_case("INS_iDEL", {"Chromosome19": "CTG"}, [{"type": "INS_iDEL", "number": 1, "length_ranges": [[1, 1], [1, 1], [1, 1]]}], ['TC', 'GT'], homozygous_only=True),
    make_case("rTRA_1", {"Chromosome19": "CTTTA"}, [{"type": "rTRA", "number": 1, "length_ranges": [[1, 1], [1, 1], [3, 3]]}], ['ATTTC'], homozygous_only=True),
    make_case("rTRA_2", {"Chromosome19": "CTTTA"}, [{"type": "rTRA", "number": 1, "length_ranges": [[1, 1], [1, 1], [3, None]]}], ['ATTTC'], homozygous_only=True),
    make_case("rTRA_3", {"Chromosome19": "CTTTA"}, [{"type": "rTRA", "number": 1, "length_ranges": [[2, 2], [2, 1], [3, None]]}], [], homozygous_only=True),
    make_case("rTRA_4", {"Chromosome19": "CTTTA"}, [{"type": "rTRA", "number": 1, "length_ranges": [[2, 2], [2, 2], [1, None]]}], ['TATCT'], homozygous_only=True),
    make_case("rTRA_5", {"Chromosome19": "CTTTA"}, [{"type": "rTRA", "number": 1, "length_ranges": [[3, 3], [2, 2], [None, None]]}], ['TACTT', 'TTACT'], homozygous_only=True),
    make_case("rTRA_6", {"Chromosome19": "CTTTA"}, [{"type": "rTRA", "number": 1, "length_ranges": [[2, 2], [2, 2], [None, None]]}], ['TATCT', 'TTCTA', 'CTATT'], homozygous_only=True),
    make_case("rTRA_7", {"ChromA": "CTTTA"}, [{"type": "rTRA", "number": 1, "length_ranges": [[2, 2], [2, 2], [None, None]]}], ['TATCT', 'TTCTA', 'CTATT'], homozygous_only=True),
    make_case("snp_A_homo", {"ChromA": "A"}, [{"type": "SNP", "number": 1}], ['C', 'T', 'G'], homozygous_only=True),
    make_case("snp_A_hetero", {"ChromA": "A"}, [{"type": "SNP", "number": 1}], ['A', 'C', 'T', 'G'], homozygous_only=False, hetero=True)
]

NO_DIS_DATA = [
    make_case("no_dis_0", {"chr21": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}, [{"type": "delINVdup", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]], "blacklist_region_type": "AluSz6"}], max_tries=2000, blacklist_regions=[EXCLUDE_BED]),
    make_case("no_dis_1", {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}, [{"type": "delINVdup", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]]}, {"type": "delINVdel", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]]}, {"type": "dupINVdup", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]]}]),
    make_case("no_dis_2", {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGAGTCAGGGAGCAAAAAAGTGTGACACTAGTCCACAGGTGAGAAACACAAATATTCAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}, [{"type": "dupINVdel", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]]}, {"type": "delINV", "number": 1, "length_ranges": [[5, 5], [5, 5]]}, {"type": "INVdel", "number": 1, "length_ranges": [[5, 5], [5, 5]]}]),
    make_case("no_dis_3", {"Chromosome19": "ACACTAGTCCACAGGTGAGAATCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}, [{"type": "dupINV", "number": 1, "length_ranges": [[5, 5], [5, 5]]}, {"type": "INVdup", "number": 1, "length_ranges": [[5, 5], [5, 5]]}]),
    make_case("no_dis_4", {"chr19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}, [{"type": "delINVdup", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]], "blacklist_region_type": "all"}], blacklist_regions=AVOID_INTERVAL),
    make_case("no_dis_5", {"Chromosome19": "CTCCGT"}, [{"type": "dupINVdup", "number": 1, "length_ranges": [[2, 2], [2, 2], [2, 2]]}], ["CTACGGAGGT"]),
    make_case("no_dis_6", {"Chromosome19": "CTCCGT"}, [{"type": "delINVdel", "number": 1, "length_ranges": [[2, 2], [2, 2], [2, 2]]}], ["GG"]),
    make_case("no_dis_7", {"Chromosome19": "CTCCGT"}, [{"type": "delINVdup", "number": 1, "length_ranges": [[2, 2], [2, 2], [2, 2]]}], ["ACGGGT"]),
    make_case("no_dis_8", {"Chromosome19": "CTCCGT"}, [{"type": "dupINVdel", "number": 1, "length_ranges": [[2, 2], [2, 2], [2, 2]]}], ["CTGGAG"]),
    make_case("no_dis_9", {"Chromosome19": "CTCCGT"}, [{"type": "delINV", "number": 1, "length_ranges": [[3, 3], [3, 3]]}], ["ACG"]),
    make_case("no_dis_10", {"Chromosome19": "CTCCGT"}, [{"type": "INVdel", "number": 1, "length_ranges": [[3, 3], [3, 3]]}], ["GAG"]),
    make_case("no_dis_11", {"Chromosome19": "CGT"}, [{"type": "DUP_INV", "number": 1, "length_ranges": [[3, 3]]}], ["ACGACG"]),
    make_case("no_dis_12", {"chr19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}, [{"type": "delINVdup", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]], "blacklist_region_type": "all"}], blacklist_regions=AVOID_INTERVAL)
]

FILTER_CHROM_DATA = [
    make_case("filter_chr_10", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA", "chr20": "CTCCGT"}, [{"type": "DEL", "number": 1, "length_ranges": [[3, 3]]}], filter_small_chr=10),
    make_case("filter_chr_4", {"chr1": "CTCCGT", "chr2": "CTCCGT", "chr3": "CTCCGT", "chrM": "C"}, [{"type": "DEL", "number": 1, "length_ranges": [[3, 3]]}], filter_small_chr=4)
]

TEST_FAIL = [
    make_case("filter_chr_50", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA", "chr20": "CTCCGT"}, [{"type": "DEL", "number": 1, "length_ranges": [[3, 3]]}], filter_small_chr=50),
    make_case("intersv_dist_2", {"Chr21": utils.generate_seq(20)}, [{"type": "DUP", "number": 3, "length_ranges": [[5, 5]]}], min_intersv_dist=10, random_seed=2),
    make_case("intersv_dist_3", {"Chr21": utils.generate_seq(10)}, [{"type": "DUP", "number": 2, "length_ranges": [[5, 5]]}], min_intersv_dist=1, random_seed=2),
    make_case("custom_overlap_1", {"chr21": "GCAGACTGAC"}, [{"type": "AB->AA", "number": 1, "length_ranges": [[6, 6], [5, 5]], "blacklist_region_type": [["L1HS"]]}], overlap_regions=[OVERLAP_BED_1]),
]

CUSTOM_SV_DATA = [
    make_case("custom_sv", {"chr21": "AGACT"}, [{"type": "A -> aB", "number": 1, "length_ranges": [[5, 5], [5, 5]]}])
]

CUSTOM_OVERLAP_DATA = [
    make_case("custom_overlap_0", {"chr21": "GCAGACTGAC"}, [{"type": "AB->AA", "number": 1, "length_ranges": [[5, 5], [5, 5]], "blacklist_region_type": [["L1HS"]]}], ["GCAGAGCAGA"], blacklist_regions=[OVERLAP_BED_1], homozygous_only=True),
]

INTSV_DISTANCE_DATA = [
    make_case("intersv_dist_0", {"Chr21": utils.generate_seq(600)}, [{"type": "DEL", "number": 20, "length_ranges": [[5, 5]]}], min_intersv_dist=10, random_seed=2),
    make_case("intersv_dist_1", {"Chr21": utils.generate_seq(1000)}, [{"type": "dDUP", "number": 20, "length_ranges": [[5, 5], [5, 50]]}], min_intersv_dist=10, random_seed=2),
]

FRAG_LEVEL_OVERLAP_DATA = [
    make_case("fo_0", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->A_A", "number": 1, "length_ranges": [[None, None], [1, 1]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_1", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->A_a", "number": 1, "length_ranges": [[None, None], [1, 1]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_2", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->_A", "number": 1, "length_ranges": [[None, None], [1, 1]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_3", {"chr21": "CTTGATGATTGATTGATGAGATTGATGATA"}, [{"type": "(A)_B->A_A", "number": 1, "length_ranges": [[None, None], [2, 2], [2, 2]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact", "overlap_region_type": ["L1"]}], overlap_regions=[OVERLAP_BED_2]),
    make_case("fo_4", {"chr21": "CTTGATGATGATTGATGATTGATTGATGAA"}, [{"type": "A_(B)->A_A", "number": 1, "length_ranges": [[2, 2], [2, 2], [None, None]], "overlap_region_length_range": [2, 2], "overlap_mode": "exact", "overlap_region_type": ["L1"]}], overlap_regions=[OVERLAP_BED_2]),
    make_case("fo_5", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->A_A", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_mode": "partial", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_6", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->A_a", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_mode": "partial", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_7", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->_A", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_mode": "partial", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_8", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "A_()->A_A", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_mode": "contained", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_9", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "A_()->A_a", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_mode": "contained", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_10", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->_A", "number": 1, "length_ranges": [[2, 2], [1, 1]], "overlap_mode": "contained", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_11", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "A(B)C->b", "number": 1, "length_ranges": [[2, 2], [None, None], [2, 2]], "overlap_region_length_range": [4, 6], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_12", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)BC->b", "number": 1, "length_ranges": [[None, None], [2, 2], [2, 2]], "overlap_region_length_range": [4, 6], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_13", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "AB(C)->b", "number": 1, "length_ranges": [[2, 2], [2, 2], [None, None]], "overlap_region_length_range": [4, 6], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_14", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)B->b", "number": 1, "length_ranges": [[None, None], [2, 2]], "overlap_region_length_range": [4, 6], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_15", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "A(B)->a", "number": 1, "length_ranges": [[2, 2], [None, None]], "overlap_region_length_range": [4, 6], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_16", {"chr21": utils.generate_seq(100)}, [{"type": "A(B)_C_D -> bb_AEc_EDC", "number": 1, "length_ranges": [[2, 3], [None, None], [10, 15], [7, 10], [15, 20], [8, 10], [10, 15]], "overlap_mode": "exact", "overlap_region_length_range": [4, 10], "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("fo_17", {"chr21": "ACCTTGA"}, [{"type": "ABC -> ACACaBC", "number": 1, "length_ranges": [[2, 2], [3, 3], [2, 2]]}], ["ACGAACGAGTCTTGA"], homozygous_only=True)
]

FRAG_LEVEL_OVERLAP_UNBOUNDED_DISP_DATA = [
    make_case("frag_disp_0", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->A_A", "number": 1, "length_ranges": [[None, None], [1, None]], "overlap_mode": "exact", "overlap_region_type": ["L1HS"]}], overlap_regions=[OVERLAP_BED_1]),
    make_case("frag_disp_1", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_ -> _A", "number": 1, "length_ranges": [[None, None], [1, None]], "overlap_mode": "exact", "overlap_region_type": ["L1", "Alu"]}], overlap_regions=[OVERLAP_BED_2]),
    make_case("frag_disp_2", {"chr21": "CTTGATGATGATTGATGATA"}, [{"type": "(A)_->A_A", "number": 1, "length_ranges": [[None, None], [1, None]], "overlap_mode": "exact", "overlap_region_type": ["Alu", "L1"]}], overlap_regions=[OVERLAP_BED_2])
]

PARTIAL_OVERLAP_DATA = [
    make_case("po_0", {"chr21": "CTCCGTAGTA"}, [{"type": "DEL", "number": 1, "length_ranges": [[2, 2]], "overlap_region_type": "L1HS", "overlap_mode": "partial"}], overlap_regions=[OVERLAP_BED_12]),
    make_case("po_1", {"chr21": "CTCCGTAGTAAGTCAGGTGAGGCAG"}, [{"type": "DEL", "number": 1, "length_ranges": [[2, 6]], "overlap_region_type": "L1HS", "overlap_mode": "partial"}, {"type": "DEL", "number": 1, "length_ranges": [[None, None]], "overlap_region_length_range": [2, 6], "overlap_region_type": "L1HS", "overlap_mode": "exact"}], overlap_regions=[OVERLAP_BED_7])
]

UNBOUNDED_DISP_DATA = [
    make_case("ub_1", {"chr21": utils.generate_seq(1000)}, [{"type": "nrTRA", "number": 1, "length_ranges": [[2, 2], [1, None]]}, {"type": "dDUP", "number": 1, "length_ranges": [[2, 2], [1, None]]}, {"type": "INV_dDUP", "number": 1, "length_ranges": [[2, 2], [1, None]]}, {"type": "dDUP_iDEL", "number": 1, "length_ranges": [[2, 2], [2, 2], [1, None]]}, {"type": "INS_iDEL", "number": 1, "length_ranges": [[2, 2], [2, 2], [1, None]]}]),
    make_case("ub_2", {"chr21": utils.generate_seq(1000)}, [{"type": "A_B_C -> a_bb_C", "number": 2, "length_ranges": [[2, 2], [2, None], [2, 2], [2, None], [2, 2]]}]),
]

REQ_SPACE_DATA = [
    make_case("req_space", {"chr21": "CTCCGT"}, [{"type": "DEL", "number": 1, "length_ranges": [[9, 9]]}])
]

ROI_PLACEMENT_FAILURE_DATA = [
    make_case("roi_fail_1", {"chr21": "CTCCGTCGTA"}, [{"type": "DEL", "number": 1, "overlap_region_length_range": [10, 10], "overlap_mode": "exact"}], overlap_regions=[OVERLAP_BED_1]),
    make_case("roi_fail_2", {"chr21": "CTCCGTCGTACTAAGTCGTA"}, [{"type": "DEL", "number": 1, "overlap_region_length_range": [10, 10], "overlap_mode": "exact"}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2]),
    make_case("roi_fail_3", {"chr19": "CTCCGTCGTACTAAGTCGTA"}, [{"type": "DEL", "number": 1, "overlap_region_length_range": [10, 10], "overlap_mode": "exact"}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2]),
    make_case("roi_fail_4", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "DEL", "number": 1, "length_ranges": [[1, 10]]}, {"type": "DEL", "number": 2, "overlap_region_length_range": [1, 10], "overlap_mode": "exact", "overlap_region_type": [["L1HS"]]}, {"type": "DEL", "number": 4, "overlap_region_length_range": [1, 10], "overlap_mode": "exact", "overlap_region_type": [["ALR/Alpha"]]}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2]),
    make_case("roi_fail_5", {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "DEL", "number": 2, "overlap_region_length_range": [1, 5], "overlap_mode": "exact", "overlap_region_type": [["L1HS"]]}, {"type": "DEL", "number": 3, "overlap_region_length_range": [1, 5], "overlap_mode": "exact", "overlap_region_type": [["ALR/Alpha"]]}], overlap_regions=[OVERLAP_BED_1, OVERLAP_BED_2]),
    make_case("roi_fail_6", {"chr21": "CCTCCGTCGTACTAAGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTATCCGTCGTACTAAGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}, [{"type": "DEL", "number": 1, "overlap_region_length_range": [2, 4], "overlap_mode": "exact", "overlap_region_type": [["Alu"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [2, 4], "overlap_mode": "exact", "overlap_region_type": [["L1"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [2, 4], "overlap_mode": "exact", "overlap_region_type": [["L2"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [2, 4], "overlap_mode": "exact", "overlap_region_type": [["SVA"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [2, 4], "overlap_mode": "exact", "overlap_region_type": [["HERVK"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [6, 8], "overlap_mode": "exact", "overlap_region_type": [["Alu"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [6, 8], "overlap_mode": "exact", "overlap_region_type": [["L1"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [6, 8], "overlap_mode": "exact", "overlap_region_type": [["L2"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [6, 8], "overlap_mode": "exact", "overlap_region_type": [["SVA"]]}, {"type": "DEL", "number": 1, "overlap_region_length_range": [6, 8], "overlap_mode": "exact", "overlap_region_type": [["HERVK"]]}], overlap_regions=[OVERLAP_BED_11])
]

BLACKLIST_REGIONS_DATA = [
    make_case("bl_1", {"chr21": "TCGATCGATCGAT"}, [{"type": "DEL", "length_ranges": [[6, 6]], "blacklist_region_type": "all", "number": 1}], ['TCTCGAT', 'TCGATCG', 'TCGCGAT'], random_seed=3, blacklist_regions=[AVOID_INTERVAL_3], homozygous_only=True, min_intersv_dist=0),
    make_case("bl_2", {"chr21": "TCGATCGATCGAT"}, [{"type": "DEL", "length_ranges": [[6, 6]], "blacklist_region_type": [['TEST1', 'TEST3']], "number": 1}], ['TCTCGAT', 'TCGATCT', 'TCGCGAT', 'TCGATCG', 'TCGATAT', 'GATCGAT']  , random_seed=3, blacklist_regions=[AVOID_INTERVAL_3], homozygous_only=True, min_intersv_dist=0),
    make_case("bl_3", {"chr21": "TCGATCGATCGAA"}, [{"type": "A_->A_A", "length_ranges": [[6, 6], [5, 5]], "blacklist_region_type": [['TEST'], ['TEST1', 'TEST2', 'TEST3'], ['TESTVCF']], "number": 1}, {"type": "DEL", "length_ranges": [[3, 3]], "number": 1}], ['TCGATCGATAGATCGA'], random_seed=3, blacklist_regions=[AVOID_INTERVAL_2, AVOID_INTERVAL_3, AVOID_INTERVAL], homozygous_only=True, min_intersv_dist=0),
    make_case("bl_4", {"chr21": "TCGATCGATCGAT"}, [{"type": "DEL", "length_ranges": [[5, 10]], "overlap_mode": "containing", "blacklist_region_type": "all", "number": 1}], ['TCGATCGA', 'TCGATCG', 'TCGATC', 'TCGAT', 'TCG', 'TCGA'], random_seed=3, overlap_regions=[AVOID_INTERVAL_2], blacklist_regions=[AVOID_INTERVAL_2], homozygous_only=True, min_intersv_dist=0),
    make_case("bl_5", {"chr21": "TCGATCGATCGAT"}, [{"type": "DEL", "length_ranges": [[5, 5]], "overlap_mode": "partial", "blacklist_region_type": "all", "number": 1}], ['TCGACGAT', 'TCGATGAT', 'CGATCGAT'], random_seed=3, overlap_regions=[AVOID_INTERVAL_3], blacklist_regions=[AVOID_INTERVAL], homozygous_only=True, min_intersv_dist=0)
]

ARM_DATA = [
    make_rt_case("arm_1", "TCGAAGCT", [{"type": "DEL", "overlap_mode": 'terminal', "length_ranges": [[4, 4]]}], ["AGCT", "TCGA"], random_seed=2),
    make_rt_case("arm_2", "TCGATCGATCGATCGA", [{"type": "DUP", "overlap_mode": 'whole-chromosome', "length_ranges": [[None, None]], "n_copies": 2}], [("TCGATCGATCGATCGA", "TCGATCGATCGATCGA", "TCGATCGATCGATCGA"), "TCGATCGATCGATCGA"], random_seed=2, hetero=True, homozygous_only=False),
    make_rt_case("arm_3", "TACGTACGGCAT", [{"type": "(AB)C ->bcb", "overlap_mode": 'terminal', "length_ranges": [[4, 4], [4, 4], [4, 4]]}], ["CGTAATGCCGTA"], random_seed=2),
    make_rt_case("arm_4", ["TCGATCGATCGATCGA", "TCGA"], [{"type": "DEL", "overlap_mode": 'whole-chromosome', "length_ranges": [[None, None]]}], [("TCGATCGATCGATCGA", "TCGA"), "TCGA", "TCGATCGATCGATCGA"], random_seed=2, hetero=True, homozygous_only=False)
]

DIVERGENCE_DUP_DATA = [
    make_case("div_dup_1", {"chr21": "TC"}, [{"type": "(A)->Aa*", "divergence_prob": 0.5, "overlap_mode": "exact", "length_ranges": [[1, 1]], "number": 1}], ["TAC", "TTC", "TGC", "TCC"], random_seed=2, overlap_regions=[INSERTION_ORDER_EXACT], homozygous_only=True, min_intersv_dist=0),
    make_case("div_dup_2", {"chr21": "TC"}, [{"type": "(A)->A*a", "divergence_prob": 0.5, "overlap_mode": "exact", "length_ranges": [[1, 1]], "number": 1}], ["TAC", "CAC", "GAC", "AAC"], random_seed=2, overlap_regions=[INSERTION_ORDER_EXACT], homozygous_only=True, min_intersv_dist=0),
    make_case("div_dup_3", {"chr21": "TC"}, [{"type": "(A)->A*a", "divergence_prob": 0.5, "overlap_mode": "exact", "length_ranges": [[1, 1]], "number": 1}], ["TC", "TAC", "CAC", "GAC", "AAC"], random_seed=2, overlap_regions=[INSERTION_ORDER_EXACT], hetero=True, homozygous_only=False, min_intersv_dist=0),
    make_case("div_dup_4", {"chr21": "TC"}, [{"type": "(A)_->A_A*", "divergence_prob": 0.5, "overlap_mode": "exact", "length_ranges": [[1, 1], [1, 1]], "number": 1}], ["TCT", "TCA", "TCC", "TCG"], random_seed=2, overlap_regions=[INSERTION_ORDER_EXACT], homozygous_only=True, min_intersv_dist=0),
    make_case("div_dup_5", {"chr21": "TCG"}, [{"type": "A->AA*", "divergence_prob": 0.5, "length_ranges": [[3, 3]], "number": 1}], ["TCG" + c0 + c1 + c2 for c0 in 'TCGA' for c1 in 'TCGA' for c2 in 'TCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0)
]

SNP_OVERLAP_DATA = [
    make_case("snp_ov_1", {"chr21": "TCG"}, [{"type": "SNP", "allow_sv_overlap": True, "number": 1}, {"import": IMPORT_DEL}], ["C", "T", "A", "G"], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("snp_ov_2", {"chr21": "TC"}, [{"type": "SNP", "allow_sv_overlap": True, "number": 1}, {"import": IMPORT_INV}], ["GC", "GT", "GG", "AA", "CA", "TA"], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("snp_ov_3", {"chr21": "TC"}, [{"type": "SNP", "allow_sv_overlap": True, "number": 1}, {"type": "DUP", "length_ranges": [[2, 2]], "number": 1}], ["TATA", "TGTG", "TTTT", "ACAC", "GCGC", "CCCC", "TA", "TG", "TT", "AC", "TC", "GC", "CC", "TCTC"], random_seed=2, homozygous_only=False, hetero=True, min_intersv_dist=0),
    make_case("snp_ov_4", {"chr21": "TC"}, [{"type": "SNP", "allow_sv_overlap": True, "number": 1}, {"type": "DUP", "length_ranges": [[2, 2]], "number": 1}], ["TATA", "TGTG", "TTTT", "ACAC", "GCGC", "CCCC"], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("snp_ov_5", {"chr21": "TC"}, [{"import": IMPORT_SNP_OVERLAP}, {"type": "DUP", "length_ranges": [[2, 2]], "number": 1}], ["TATA", "TCTC"], random_seed=2, min_intersv_dist=0),
    make_case("snp_ov_6", {"chr21": "TC"}, [{"type": "SNP", "allow_sv_overlap": True, "number": 1}, {"type": "dDUP", "length_ranges": [[1, 1], [1, 1]], "number": 1}], ['T' + s + 'T' for s in 'TGA'] + ['C' + s + 'C' for s in 'CGA'] + [s + 'T' + s for s in 'TGA'] + [s + 'C' + s for s in 'CGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("snp_ov_7", {"chr21": "TC"}, [{"type": "SNP", "allow_sv_overlap": True, "number": 1}, {"type": "INV", "length_ranges": [[2, 2]], "number": 1}], ['G' + s for s in 'TCG'] + [s + 'A' for s in 'TCA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("snp_ov_8", {"chr21": "TC"}, [{"type": "SNP", "allow_sv_overlap": True, "number": 1}, {"type": "A_->A_A+", "n_copies": 3, "length_ranges": [[1, 1], [1, 1]], "number": 1}], ['T' + s + 'TTT' for s in 'TGA'] + [s + 'C' + s + s + s for s in 'CGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("snp_ov_9", {"chr21": "TCG"}, [{"type": "SNP", "length_ranges": [[1, 1]], "allow_sv_overlap": True, "number": 1}, {"type": "A_->_A+", "length_ranges": [[2, 2], [1, 1]], "number": 1, "n_copies": 3}], [template.replace(target_char, mut) for template in ["GTCTCTC", "TCG"] for target_char in set(template) for mut in "TCGA"], random_seed=2, homozygous_only=False, allow_hap_overlap=True, hetero=True, min_intersv_dist=0),
    make_case("snp_ov_10", {"chr21": "TCGA"}, [{"type": "INDEL", "length_ranges": [[2, 2]], "allow_sv_overlap": True, "number": 1}, {"type": "A_->_A", "length_ranges": [[2, 2], [2, 2]], "number": 1}], ['TC', 'AT', 'TC', 'GA'] + ['GA' + s1 + s2 + 'TC' for s1 in 'TCGA' for s2 in 'TCGA'] + [s1 + s2 + 'GATC' for s1 in 'TCGA' for s2 in 'TCGA'] + ['G' + s1 + s2 + 'ATC' for s1 in 'TCGA' for s2 in 'TCGA'] + ['GAT' + s1 + s2 + 'C' for s1 in 'TCGA' for s2 in 'TCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("snp_ov_11", {"chr21": "TCGA"}, [{"type": "INDEL", "length_ranges": [[1, 1]], "allow_sv_overlap": True, "number": 1}, {"import": IMPORT_DEL}, {"type": "INV", "length_ranges": [[2, 2]], "number": 1}], ['TC' + s for s in 'TCGA'] + [s + 'TC' for s in 'TCGA'] + ['T' + s + 'C' for s in 'TCGA'] + ['C', 'T', 'TC'], random_seed=2, homozygous_only=True, min_intersv_dist=0)
]

INDEL_OVERLAP_DATA = [
    make_case("indel_ov_1", {"chr21": "TC"}, [{"type": "(A)->Aa", "overlap_mode": "exact", "length_ranges": [[1, 1]], "number": 1}, {"type": "A->aA", "length_ranges": [[1, 1]], "number": 1}], ["TAGC"], random_seed=2, overlap_regions=[INSERTION_ORDER_EXACT], homozygous_only=True, min_intersv_dist=0),
    make_case("indel_ov_2", {"chr21": "TC"}, [{"type": "A->aA", "length_ranges": [[1, 1]], "number": 1}, {"type": "(A)->Aa", "overlap_mode": "exact", "length_ranges": [[1, 1]], "number": 1}], ["TAGC"], random_seed=2, overlap_regions=[INSERTION_ORDER_EXACT], homozygous_only=True, min_intersv_dist=0),
    make_case("indel_ov_3", {"chr21": "TCGA"}, [{"type": "INDEL", "length_ranges": [[1, 1]], "allow_sv_overlap": True, "number": 1}, {"import": IMPORT_DEL}], ['G', 'A', 'GA'] + [s + 'GA' for s in 'TCGA'] + ['GA' + s for s in 'TCGA'] + ['G' + s + 'A' for s in 'TCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("indel_ov_4", {"chr21": "TC"}, [{"type": "INDEL", "length_ranges": [[1, 1]], "allow_sv_overlap": True, "number": 1}, {"type": "DUP", "length_ranges": [[2, 2]], "number": 1}], ['TT', 'CC'] + [s + 'TCTC' for s in 'TCGA'] + ['T' + s + 'CT' + s + 'C' for s in 'TCGA'] + ['TCTC' + s for s in 'TCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("indel_ov_5", {"chr21": "TC"}, [{"type": "DUP", "length_ranges": [[2, 2]], "number": 1}, {"type": "INDEL", "length_ranges": [[1, 1]], "allow_sv_overlap": True, "number": 1}], ['TT', 'CC'] + [s + 'TCTC' for s in 'TCGA'] + ['T' + s + 'CT' + s + 'C' for s in 'TCGA'] + ['TCTC' + s for s in 'TCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("indel_ov_6", {"chr21": "TC"}, [{"type": "INDEL", "length_ranges": [[1, 1]], "allow_sv_overlap": True, "number": 1}, {"type": "DUP", "n_copies": 3, "length_ranges": [[2, 2]], "number": 1}], ['TTTT', 'CCCC'] + [s + 'TCTCTCTC' for s in 'TCGA'] + ['T' + s + 'CT' + s + 'CT' + s + 'CT' + s + 'C' for s in 'TCGA'] + ['TCTCTCTC' + s for s in 'TCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("indel_ov_7", {"chr21": "TCG"}, [{"type": "INDEL", "length_ranges": [[1, 1]], "allow_sv_overlap": True, "number": 1}, {"type": "A_->_A", "length_ranges": [[2, 2], [1, 1]], "number": 1}], ['GC', 'GT', 'TC'] + ['G' + s + 'TC' for s in 'TCGA'] + [s + 'GTC' for s in 'TCGA'] + ['GT' + s + 'C' for s in 'TCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("indel_ov_8", {"chr21": "TCGA"}, [{"type": "INDEL", "length_ranges": [[2, 2]], "allow_sv_overlap": True, "number": 1}, {"type": "A_->_A", "length_ranges": [[2, 2], [2, 2]], "number": 1}], ['TC', 'AT', 'TC', 'GA'] + ['GA' + s1 + s2 + 'TC' for s1 in 'TCGA' for s2 in 'TCGA'] + [s1 + s2 + 'GATC' for s1 in 'TCGA' for s2 in 'TCGA'] + ['G' + s1 + s2 + 'ATC' for s1 in 'TCGA' for s2 in 'TCGA'] + ['GAT' + s1 + s2 + 'C' for s1 in 'TCGA' for s2 in 'TCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("indel_ov_9", {"chr21": "TCGA"}, [{"type": "INDEL", "length_ranges": [[1, 1]], "allow_sv_overlap": True, "number": 1}, {"import": IMPORT_DEL}, {"type": "INV", "length_ranges": [[2, 2]], "number": 1}], ['TC' + s for s in 'TCGA'] + [s + 'TC' for s in 'TCGA'] + ['T' + s + 'C' for s in 'TCGA'] + ['C', 'T', 'TC'], random_seed=2, homozygous_only=True, min_intersv_dist=0)
]

SIMPLE_TR_DATA = [
    make_case("tr_1", {"chrA": "TCGTCGCGGATATAT"}, [{"type": "trCON", "repeat_count_change_range": [3, 3], "number": 1}, {"type": "trCON", "repeat_count_change_range": [2, 2], "number": 1}], ["TCGTG", "CGG"], random_seed=2, overlap_regions=[MOTIF_BED], homozygous_only=True, min_intersv_dist=0),
    make_case("tr_2", {"chrA": "TCGTCGCGGATATAT"}, [{"type": "trEXP", "repeat_count_change_range": [3, 3], "number": 1}, {"type": "trEXP", "repeat_count_change_range": [2, 2], "number": 1}], ["TCGTCGTCGTCGCGGATATATATATAT", "TCGTCGCGCGCGGATATATATATAT", "TCGTCGCGCGCGCGGATATATATAT", "TCGTCGTCGTCGTCGCGGATATATATAT"], random_seed=2, overlap_regions=[MOTIF_BED], homozygous_only=True, min_intersv_dist=0),
    make_case("tr_3", {"chrA": "TCGTCGCGGATATAT"}, [{"type": "trCON", "repeat_count_change_range": [3, 3], "number": 1}, {"type": "trEXP", "repeat_count_change_range": [2, 2], "number": 1}], ["TCGTCGTCGTCGCGG", "TCGTCGCGCGCGG"], random_seed=2, overlap_regions=[MOTIF_BED], homozygous_only=True, min_intersv_dist=0)
]

SMALL_CHR_FILTER_DATA = [
    make_rt_case("small_chr_1", ["TCG", "TCGATCGATCGA"], [{"type": "A->", "length_ranges": [[2, 2]]}], [("TCG", "TCGATCGATCGA"[:i] + ("TCGATCGATCGA"[i+2:] if i < 11 else '')) for i in range(11)], random_seed=2, filter_small_chr=4)
]

MULTIPLE_OV_BLCK_FILE_DATA = [
    make_case("mob_1", {"chr21": "TCGATCGATCGATCGA"}, [{"type": "DEL", "number": 1, "overlap_mode": "contained", "overlap_region_type": [["all"], ["Alu"]], "length_ranges": [[4, 4]]}], ["TCGATCGATCGA", "TCGATCGATCGA", "TCGATCGATCGA"], random_seed=2, overlap_regions=[OVERLAP_VCF, OVERLAP_BED_2], homozygous_only=True, min_intersv_dist=0),
    
    make_case("mob_2", {"chr21": "TCGATCGATCGATCGA"}, [{"type": "DEL", "number": 1, "overlap_mode": "exact", "overlap_region_length_range": [None, 1], "overlap_region_type": [["all"], ["Alu"]], "length_ranges": [[None, None]]}], ["TGATCGATCGATCGA"], random_seed=2, overlap_regions=[OVERLAP_VCF, OVERLAP_BED_2], homozygous_only=True, min_intersv_dist=0),
    
    make_case("mob_3", {"chr21": "TCGATCGATCGATCGA"}, [{"type": "DEL", "number": 1, "overlap_mode": "contained", "overlap_region_type": [["TEST2"], ["Alu"]], "length_ranges": [[3, 3]]}], ["TCGATCCGATCGA", "TCGATCGGATCGA", "TCGATCGAATCGA", "TCGATCGATTCGA"], random_seed=2, overlap_regions=[OVERLAP_VCF, OVERLAP_BED_2], homozygous_only=True, min_intersv_dist=0),
    
    make_case("mob_4", {"chr21": "TCGATCGATCGATCGA"}, [{"type": "DEL", "number": 1, "overlap_mode": "exact", "overlap_region_type": [["TEST", "TEST2"], ["Alu", "L1"]], "length_ranges": [[None, None]]}], ["TGATCGATCGATCGA", "TCGCGATCGATCGA", "TCGATCTCGA", "TCGATCGATCGATA"], random_seed=2, overlap_regions=[OVERLAP_VCF, OVERLAP_BED_2], homozygous_only=True, min_intersv_dist=0),
    
    make_case("mob_5", {"chr21": "TCGATCGATCGATCGA"}, [{"type": "DEL", "number": 1, "blacklist_region_type": [["TEST", "TEST2"], ["Alu", "L1"]], "length_ranges": [[3, 3]]}], ["ATCGATCGATCGA", "TCGATCGATCGAT", "TCCGATCGATCGA", "TCGATCGATCGAA", "TCGGATCGATCGA"], random_seed=2, blacklist_regions=[OVERLAP_VCF, OVERLAP_BED_2], homozygous_only=True, min_intersv_dist=0)
]

SIMPLE_TEST_DATA = [
    make_rt_case("st_1", "TCT", [{"type": "AB->A+B+", "n_copies": [3, 5], "length_ranges": [[2, 2], [1, 1]]}], ["TCTCTCTTTTT"]),
    make_rt_case("st_2", "TCTAGTCCTGGTAAT", [{"type": "_AB_C_->c_AB_AC_B", "length_ranges": [[2, 2], [2, 2], [1, 1], [1, 1], [3, 3], [3, None]]}], ["AGGTCTAGTTACCTGGTGAAT", "AGGTCTAGTTACCTGGTAGAT", "AGGTCTAGTTACCTGGTAAGT", "AGGTCTAGTTACCTGGTAATG", "TCAGCTAGTCAGCTGGTATAT", "TCAGCTAGTCAGCTGGTAATT", "TCCCATAGTCCGTTGGTAACT", "TCCCATAGTCCGTTGGTAATC", "TCTACCAGTCCTTCGGTAATC"]),
    make_rt_case("st_3", "TC", [{"type": "A->A+", "n_copies": [3], "length_ranges": [[2, 2]]}], ["TCTCTC"]),
    make_rt_case("st_4", "TCA", [{"type": "A->", "length_ranges": [[2, 2]]}], ["T", "A"]),
    make_rt_case("st_5", "TCA", [{"type": "A  ->  ", "length_ranges": [[2, 2]]}], ["T", "A"]),
    make_rt_case("st_6", "TCA", [{"type": "mCNV", "n_copies": [3], "n_copiesB": [0], "length_ranges": [[2, 2]]}], ["TCTCTCA", "T", "A", "TCACACA"], homozygous_only=False),
    make_rt_case("st_7", "TCA", [{"type": "mCNV", "n_copies": [[2, 3]], "n_copiesB": [0], "length_ranges": [[2, 2]]}], ["TCTCTCA", "TCTCA", "T", "A", "TCACACA", "TCACA"], homozygous_only=False),
    make_rt_case("st_8", "TCA", [{"type": "mCNV", "n_copies": [[2, 3]], "n_copiesB": [[2, 2]], "length_ranges": [[2, 2]]}], ["TCTCTCA", "TCTCA", "TCACACA", "TCACA"], homozygous_only=False),
    make_rt_case("st_9", "TC", [{"type": "A->A+", "n_copies": [[2, 3]], "length_ranges": [[2, 2]]}], ["TCTCTC", "TCTC"]),
    make_rt_case("st_10", "TCGA", [{"type": "AB->ABA+", "n_copies": [3], "length_ranges": [[2, 2], [2, 2]]}], ["TCGATCTCTC"]),
    make_rt_case("st_11", "TCGA", [{"type": "AB->A+BA", "n_copies": [3], "length_ranges": [[2, 2], [2, 2]]}], ["TCTCTCGATC"]),
    make_rt_case("st_12", "TCGA", [{"type": "AB->A+BAB+", "n_copies": [3, [1, 2]], "length_ranges": [[2, 2], [2, 2]]}], ["TCTCTCGATCGA", "TCTCTCGATCGAGA"]),
    make_rt_case("st_13", "TC", [{"type": "INS", "length_ranges": [[1, 1]]}], ["TC" + ins for ins in 'TCGA'] + [ins + "TC" for ins in 'TCGA'] + ["T" + ins + "C" for ins in 'TCGA']),
    make_rt_case("st_14", "TC", [{"type": "A->B", "novel_insertions": "tests/inputs/novel_insertions.bed", "length_ranges": [[2, 2], [None, None]]}], ["GAG"]),
    make_rt_case("st_15", "TC", [{"type": "A->BB", "novel_insertions": "tests/inputs/novel_insertions.bed", "length_ranges": [[2, 2], [None, None]]}], ["GAGGAG"]),
    make_rt_case("st_16", "TC", [{"type": "A->BAB", "novel_insertions": "tests/inputs/novel_insertions.bed", "length_ranges": [[2, 2], [3, 3]]}], ["GAGTCGAG"]),
    make_rt_case("st_17", "TCGA", [{"type": "A_->B_A", "novel_insertions": "tests/inputs/novel_insertions.bed", "length_ranges": [[2, 2], [2, 2], [3, 3]]}], ["GAGGATC"]),
    make_rt_case("st_18", "TCGA", [{"type": "DEL", "length_ranges": [[2, 2]]}], ["TC", "TA", "GA"]),
    make_rt_case("st_19", ["TCGA", "CTG"], [{"type": "DEL", "length_ranges": [[2, 2]]}], [("TC", "CTG"), ("TA", "CTG"), ("GA", "CTG"), ("TCGA", "C"), ("TCGA", "G")]),
    make_rt_case("st_20", "TCGA", [{"type": "DEL", "length_ranges": [[3, 3]]}], ["T", "A"]),
    make_rt_case("st_21", "TC", [{"type": "SNP", "length_ranges": [[1, 1]]}], ["TT", "TG", "TA", "AC", "GC", "CC"]),
    make_rt_case("st_22", "T", [{"type": "SNP", "length_ranges": [[1, 1]]}], ["C", "G", "A"]),
    make_rt_case("st_23", "TCGA", [{"type": "INV", "length_ranges": [[3, 4]]}], ["TCGA", "TTCG", "CGAA"]),
    make_rt_case("st_24", "TCGA", [{"type": "DUP", "length_ranges": [[3, 4]]}], ["TCGTCGA", "TCGACGA", "TCGATCGA"]),
    make_rt_case("st_25", "TCGA", [{"type": "nrTRA", "length_ranges": [[2, 2], [1, None]]}], ["CGTA", "TACG", "GATC", "GTCA", 'TGAC']),
    make_rt_case("st_26", "TCGA", [{"type": "rTRA", "length_ranges": [[2, 2], [2, 2], [None, None]]}], ["GATC"]),
    make_rt_case("st_27", "TCGA", [{"type": "rTRA", "length_ranges": [[2, 2], [1, 1], [None, None]]}], ["AGTC", "CGTA", "TACG", "GACT", "TGAC", 'GTCA']),
    make_rt_case("st_28", "TCGA", [{"type": "rTRA", "length_ranges": [[2, 2], ["A", "A"], [None, None]]}], ["GATC"]),
    make_rt_case("st_29", "TCGA", [{"type": "rTRA", "length_ranges": [["B", "B"], [2, 2], [None, None]]}], ["GATC"]),
    make_rt_case("st_30", "TCGA", [{"type": "rTRA", "length_ranges": [[3, 3], ["A-2", "A-2"], [None, None]]}], ["ATCG", "CGAT"]),
    make_rt_case("st_31", "TCGA", [{"type": "rTRA", "length_ranges": [[1, 1], [1, 1], [2, None]]}], ["ACGT"]),
    make_rt_case("st_32", "TCGA", [{"type": "dupINVdup", "length_ranges": [[1, 1], [2, 2], [1, 1]]}], ["TTCGAA"]),
    make_rt_case("st_33", "TCGA", [{"type": "delINVdel", "length_ranges": [[1, 1], [2, 2], ["A", "A"]]}], ["CG"]),
    make_rt_case("st_34", "TCGA", [{"type": "delINVdup", "length_ranges": [[1, 1], [2, 2], [1, 1]]}], ["TCGA"]),
    make_rt_case("st_35", "TCGAG", [{"type": "delINVdup", "length_ranges": [[1, 1], [3, 3], [1, 1]]}], ["CTCGG"]),
    make_rt_case("st_36", "TCGA", [{"type": "dupINVdel", "length_ranges": [[1, 1], [2, 2], [1, 1]]}], ["TCGA"]),
    make_rt_case("st_37", "TCGAG", [{"type": "dupINVdel", "length_ranges": [[1, 1], [3, 3], [1, 1]]}], ["TTCGA"]),
    make_rt_case("st_38", "TCGA", [{"type": "delINV", "length_ranges": [[2, 2], [2, 2]]}], ["TC"]),
    make_rt_case("st_39", "TCGA", [{"type": "INVdel", "length_ranges": [[2, 2], [2, 2]]}], ["GA"]),
    make_rt_case("st_40", "TCGA", [{"type": "INS_iDEL", "length_ranges": [[2, 2], [1, 1], [1, None]]}], ["GTC", "GAC"]),
    make_rt_case("st_41", "TCGA", [{"type": "DUP_INV", "length_ranges": [[2, 2]]}], ["GAGAGA", "TCTCTC", "TCGCGA"]),
    make_rt_case("st_42", "TCGA", [{"type": "INV_DUP", "length_ranges": [[3, 3]]}], ["TCGCGAA", "TCGATCG", "CGATCGA", "TTCGCGA"]),
    make_rt_case("st_43", "TCGA", [{"type": "dupINV", "length_ranges": [[2, 2], [2, 2]]}], ["TCTCGA"]),
    make_rt_case("st_44", "TCGA", [{"type": "INVdup", "length_ranges": [[2, 2], [2, 2]]}], ["TCGAGA"]),
    make_rt_case("st_45", ["TCGA", "CAT"], [{"type": "dDUP", "length_ranges": [[3, 3], [1, None]]}], [("TCGATCG", "CAT"), ("CGATCGA", "CAT")]),
    make_rt_case("st_46", ["TCGA", "CAT"], [{"type": "dDUP", "length_ranges": [[3, 3], [None, None]], "interchromosomal": True}], [("TCGA"[:i] + "CAT" + "TCGA"[i:], "CAT") for i in range(0, len("TCGA") + 1)] + [("TCGA", "CAT"[:j] + s + "CAT"[j:]) for j in range(0, len("CAT") + 1) for s in ("TCG", "CGA")]),
    make_rt_case("st_47", ["TCGA", "CAT"], [{"type": "INV_dDUP", "length_ranges": [[3, 3], [1, None]]}], [('TCGACGA', 'CAT'), ('TCGTCGA', 'CAT')]),
    make_rt_case("st_48", ["TCGA", "CAT"], [{"type": "INV_dDUP", "length_ranges": [[3, 3], [None, None]]}], [('TCGATCG', 'CAT'), ('TTCGCGA', 'CAT'), ('TCGTCGA', 'CAT'), ('TCGACGA', 'CAT'), ('TCGCGAA', 'CAT'), ('CGATCGA', 'CAT'), ('TCGA', 'CATATG'), ('TCGA', 'ATGCAT')], homozygous_only=False, hetero=True),
    make_rt_case("st_49", "TCGA", [{"type": "INV_nrTRA", "length_ranges": [[3, 3], [1, 1]]}], ["ACGA", "TCGT"]),
    make_rt_case("st_50", ["GGTT", "CA"], [{"type": "INV_rTRA", "length_ranges": [[3, 3], [1, 1], [None, None]], "interchromosomal_period": 0}], [('GT', 'ACCA'), ('GT', 'CAAC'), ('GG', 'AACA'), ('TT', 'CACC')]),
    make_rt_case("st_51", "T", [{"type": "A->AA*", "length_ranges": [[1, 1]], "divergence_prob": [0.5]}], ["TT", "TC", "TG", "TA"]),
    make_rt_case("st_52", ["GGCCTT", "CA"], [{"type": "ABC->AC", "length_ranges": [[2, 2], [2, 2], [2, 2]]}], [("GGTT", "CA")]),
    make_rt_case("st_53", ["GGTT", "CA"], [{"type": "rTRA", "interchromosomal_period": 0, "length_ranges": [[1, 4], [1, 2], [None, None]]}], [("CGTT", "GA"), ("GGCT", "TA"), ("CTT", "GGA"), ("GGTA", "CT"), ("GGCA", "TT"), ('CA', 'GGTT'), ('CAT', 'GGT'), ('GGA', 'CTT'), ('GGTCA', 'T'), ('GGCAT', 'T'), ('GCATT', 'G'), ('GA', 'CGTT'), ('CAGTT', 'G'), ('A', 'CGGTT'), ('GATT', 'CG'), ('AGTT', 'CG'), ('GGAT', 'CT'), ('CT', 'GGTA'), ('AT', 'CGGT'), ('CT', 'GGTA'), ('GCA', 'GTT'), ('CATT', 'GG'), ('GCTT', 'GA'), ('GC', 'GTTA'), ('C', 'GGTTA'), ('GCAT', 'GT'), ('GGTC', 'TA'), ('GGC', 'TTA'), ('GCT', 'GTA'), ('ATT', 'CGG'), ('GAT', 'CGT')]),
    make_rt_case("st_54", ["GGCCTTG", "CA"], [{"type": "A_B_C->A__C", "length_ranges": [[1, 1], [1, 1], [1, 1], [2, 2], [1, 1]]}, {"type": "rTRA", "interchromosomal_period": 0, "length_ranges": [[1, 1], [2, 2], [None, None]]}], [('GGCAG', 'CTT'), ('GGCCG', 'TTA'), ('CAGCTTG', 'G'), ('GGCTG', 'CTA'), ('GGCCATG', 'T'), ('GGCTTCA', 'G'), ('GGCTCAG', 'T'), ('GCACTTG', 'G'), ('GGCATTG', 'C'), ('GGATG', 'CCT')]),
    make_rt_case("st_55", "TCGA", [{"type": "ABCD->cbda", "length_ranges": [[1, 1], [1, 1], [1, 1], [1, 1]]}], ["CGTA"]),
    make_rt_case("st_56", "GGACCT", [{"type": "ABC__->A_B_C", "length_ranges": [[1, 1], [2, 2], [1, 1], [1, 1], [1, 1]]}], ['GCGATC']),
    make_rt_case("st_57", "GGACCT", [{"type": "__ABC->A_B_C", "length_ranges": [[1, 1], [1, 1], [1, 1], [2, 2], [1, 1]]}], ['AGCCGT']),
    make_rt_case("st_58", "GGACCT", [{"type": "_ABC_->A_B_C", "length_ranges": [[1, 1], [1, 1], [2, 2], [1, 1], [1, 1]]}], ['GGACTC']),
    make_case("st_59", {"chr21": "TCGA"}, [{"import": IMPORT_SNP}], ['ACTA', 'ACGC'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_case("st_60", {"chr1": "TCGATCGA"}, [{"import": IMPORT_TEST}], ['TCGAGACGTTCG', 'TCGATCCGA'], random_seed=2, homozygous_only=True, min_intersv_dist=0),
    make_rt_case("st_61", "TCGA", [{"type": "ABCD->ABCDBABDC", "length_ranges": [[1, 1], [1, 1], [1, 1], [1, 1]]}], ['TCGACTCAG'], random_seed=2),
]

INTERCHROM_PERIOD_DATA = [
    make_rt_case("ic_1", ["TCGA", "AG"], [{"type": "A_B_C_->A_B_C_ABC", "interchromosomal_period": 1, "length_ranges": [[2, 2], [None, None], [2, 2], [None, None], [2, 2], [None, None]]}], [("TCGA", "AG"[:insertion_idx] + A + "AG" + C + "AG"[insertion_idx:]) for A in ['TC', 'GA'] for C in ['TC', 'GA'] for insertion_idx in [0, 2] if A != C], random_seed=2),
    make_rt_case("ic_2", ["TC", "AC", "TG"], [{"type": "A_B_C_->A_B_C_ABC", "interchromosomal_period": 2, "length_ranges": [[2, 2], [None, None], [2, 2], [None, None], [2, 2], [None, None]]}], 
                 [("TC", "AC"[:insertion_idx] + A + B + C + "AC"[insertion_idx:], "TG") for A in ['AC'] for C in ['TC', 'TG'] for B in ['TC', 'TG'] for insertion_idx in [0, 2] if C != B] +
                 [("TC"[:insertion_idx] + A + B + C + "TC"[insertion_idx:], "AC", "TG") for A in ['TC'] for C in ['AC', 'TG'] for B in ['AC', 'TG'] for insertion_idx in [0, 2] if C != B] +
                 [("TC", "AC", "TG"[:insertion_idx] + A + B + C + "TG"[insertion_idx:]) for A in ['TG'] for C in ['TC', 'AC'] for B in ['TC', 'AC'] for insertion_idx in [0, 2] if C != B], random_seed=2)
]


# ==========================================
# TESTS
# ==========================================

@pytest.mark.parametrize("case", SNPS_DATA + SIMPLE_DELS_DATA + SIMPLE_DUPS_DATA + SIMPLE_INVS_DATA + MULTI_INS_DATA + SIMPLE_TEST_DATA, ids=lambda c: c.id)
def test_simple_variants(case, tmp_path):
    if case.id == "snp_multi":
        frag1, frag2, sim = run_single_sv_case(case, tmp_path)
        assert frag1 not in case.ref_seqs["chr21"] or frag2 not in case.ref_seqs["chr21"]
        assert len(frag1) == len(frag2) == 10
    else:
        run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", DISPERSION_DATA, ids=lambda c: c.id)
def test_dispersions(case, tmp_path):
    run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", SIMPLE_INSS_DATA, ids=lambda c: c.id)
def test_simple_insertions(case, tmp_path):
    frag1, frag2, sim = run_single_sv_case(case, tmp_path)
    hap_bools = [len(frag) == 7 and (frag[0] == 'C' or frag[-1] == 'A') for frag in [frag1, frag2]]
    assert any(hap_bools)


@pytest.mark.parametrize("idx, case", enumerate(OVERLAP_SIMPLE_DATA))
def test_overlap_placement_simple(idx, case, tmp_path):
    frag1, frag2, _ = run_single_sv_case(case, tmp_path)
    if idx == 0:
        assert 'CTGTCGTA' in [frag1, frag2]
    elif idx == 1:
        assert 'CCCC' in frag1 or 'CCCC' in frag2
    elif idx == 2:
        assert 'AA' not in frag1 or 'AA' not in frag2


@pytest.mark.parametrize("idx, case", enumerate(COMPLEX_OVERLAP_DATA))
def test_overlap_placement_complex(idx, case, tmp_path):
    if idx <= 11:
        frag1, frag2, svs = run_single_sv_case(case, tmp_path)
        if idx == 0:
            assert frag1 in ['CTGATGA', 'CGATGAT'] or frag2 in ['CTGATGA', 'CGATGAT']
        elif idx == 1:
            assert frag1[:4] in ['CTGA', 'ACTG'] or frag2[:4] in ['CTGA', 'ACTG']
            assert frag1[-7:] in ['TCATGGA', 'ATGGATC'] or frag2[-7:] in ['TCATGGA', 'ATGGATC']
        elif idx == 2:
            for sv in svs:
                assert sv_source_segments(sv) in [[(10 + i * 3, 13 + i * 3), (13 + i * 3, 16 + i * 3), (16 + i * 3, 19 + i * 3)] for i in range(3)]
        elif idx == 3:
            for sv in svs:
                grammar = sv.config_descr.split('type: ')[-1]
                if grammar == '(A)B->b':
                    assert sv_source_segments(sv) in [[(15 + i * 3, 18 + i * 3), (18 + i * 3, 21 + i * 3)] for i in range(3)]
                elif grammar == '(A)B->a':
                    assert sv_source_segments(sv) in [[(8 + i * 2, 10 + i * 2), (10 + i * 2, 12 + i * 2)] for i in range(3)]
        elif idx == 4:
            possible_shifts = [(10 + i * 3, 13 + i * 3, 16 + i * 3) for i in range(3)]
            correct_shift = False
            (frag_a, frag_b, frag_c) = tuple(seg[0] for seg in sv_source_segments(svs[0]))
            for (ivl_a, ivl_b, ivl_c) in possible_shifts:
                if all(np.abs(frag - ivl) < 3 for (frag, ivl) in zip((ivl_a, ivl_b, ivl_c), (frag_a, frag_b, frag_c))):
                    correct_shift = True
            assert correct_shift
        elif idx in [5, 6]:
            dispersion_target = [op.target_region for op in svs[0].operations if not op.transform.is_in_place][0]
            assert 13 <= dispersion_target.start <= 16
        elif idx == 7:
            placement = [op.placement for op in svs[0].operations if not op.transform.is_in_place][0]
            breakend_pos = sorted([locus.pos for locus in placement.values()])
            assert (13 <= breakend_pos[0] <= 14) or (14 <= breakend_pos[-1] <= 15)
        elif idx == 8:
            inv_ddup = [sv for sv in svs if "A_()->A_a" in sv.info.get('GRAMMAR', '')][0]
            dispersion_target = [op.target_region for op in inv_ddup.operations if not op.transform.is_in_place][0]
            assert (10 <= dispersion_target.start <= 13) or (16 <= dispersion_target.start <= 20)
        elif idx == 9:
            grammar = svs[0].info.get('GRAMMAR', '')
            segs = sv_source_segments([sv for sv in svs if "(ABC)->b" in sv.info.get('GRAMMAR', grammar)][0])
            assert segs == [(12, 14), (14, 17), (17, 18)] or segs == [(13, 15), (15, 18), (18, 19)]
        elif idx == 10:
            segs = sv_source_segments(svs[0])
            assert segs in [[(12, 15), (15, 18)], [(13, 16), (16, 19)]]
        elif idx == 11:
            frag_bounds = sv_source_segments(svs[0])
            assert is_overlapping([(12, 19)], (frag_bounds[0][0], frag_bounds[-1][1]))
    else:
        run_seed_search(case, tmp_path)


@pytest.mark.parametrize("idx, case", enumerate(FRAG_LEVEL_OVERLAP_DATA))
def test_frag_level_overlap(idx, case, tmp_path):
    if idx <= 10:
        frag1, frag2, svs = run_single_sv_case(case, tmp_path)
        disp_ev = utils.Region(chrom='chr21', start=svs[0].placement[1].pos, end=svs[0].placement[2].pos)
    elif case.id == "fo_17":
        run_seed_search(case, tmp_path)
    else:
        run_single_sv_case(case, tmp_path)


@pytest.mark.parametrize("idx, case", enumerate(PARTIAL_OVERLAP_DATA))
def test_partial_overlap_placement(idx, case, tmp_path):
    frag1, frag2, svs = run_single_sv_case(case, tmp_path)
    svs.sort(key=lambda x: get_span(x).start)
    if idx == 0:
        assert any(p not in frag1 for p in ['CT', 'CC', 'CG', 'GT']) or any(p not in frag2 for p in ['CT', 'CC', 'CG', 'GT'])
    if idx == 1:
        case_a = (get_span(svs[0]).start, get_span(svs[0]).end) == (2, 4) and is_overlapping([(15, 20)], (get_span(svs[1]).start, get_span(svs[1]).end))
        case_b = (get_span(svs[1]).start, get_span(svs[1]).end) == (15, 20) and is_overlapping([(2, 4)], (get_span(svs[0]).start, get_span(svs[0]).end))
        assert (case_a or case_b)


@pytest.mark.parametrize("case", REQ_SPACE_DATA + ROI_PLACEMENT_FAILURE_DATA, ids=lambda c: c.id)
def test_placement_failures(case, tmp_path):
    with pytest.raises(Exception):
        run_single_sv_case(case, tmp_path)


@pytest.mark.parametrize("case", NO_DIS_DATA + FILTER_CHROM_DATA + CUSTOM_SV_DATA + CUSTOM_OVERLAP_DATA + INTSV_DISTANCE_DATA + FRAG_LEVEL_OVERLAP_UNBOUNDED_DISP_DATA, ids=lambda c: c.id)
def test_general_simulator_cases(case, tmp_path):
    if case.expected_outputs:
        run_seed_search(case, tmp_path)
    else:
        run_single_sv_case(case, tmp_path)


@pytest.mark.parametrize("case", TEST_FAIL)
def run_expected_failure_case(case: SVTestCase, tmp_path):
    """
    Runs a configuration that is expected to fail (e.g., all chromosomes filtered out).
    The test passes if and only if SVSimulator fails and raises an Exception.
    """
    with pytest.raises(Exception):
        run_single_sv_case(case, tmp_path)


@pytest.mark.parametrize("case", UNBOUNDED_DISP_DATA, ids=lambda c: c.id)
def test_unbounded_dispersion(case, tmp_path):
    run_single_sv_case(case, tmp_path)


@pytest.mark.parametrize("case", BLACKLIST_REGIONS_DATA, ids=lambda c: c.id)
def test_blacklist_regions(case, tmp_path):
    run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", INTERCHROM_PERIOD_DATA, ids=lambda c: c.id)
def test_interchromosomal_period(case, tmp_path):
    run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", ARM_DATA, ids=lambda c: c.id)
def test_simple_gain_loss(case, tmp_path):
    run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", DIVERGENCE_DUP_DATA, ids=lambda c: c.id)
def test_divergence(case, tmp_path):
    run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", SNP_OVERLAP_DATA, ids=lambda c: c.id)
def test_snp_overlap(case, tmp_path):
    run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", INDEL_OVERLAP_DATA, ids=lambda c: c.id)
def test_indel_overlap(case, tmp_path):
    run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", SIMPLE_TR_DATA, ids=lambda c: c.id)
def test_tr_operations(case, tmp_path):
    run_seed_search(case, tmp_path, allow_fail=True)


@pytest.mark.parametrize("case", SMALL_CHR_FILTER_DATA, ids=lambda c: c.id)
def test_small_chr_filter(case, tmp_path):
    run_seed_search(case, tmp_path)


@pytest.mark.parametrize("case", MULTIPLE_OV_BLCK_FILE_DATA, ids=lambda c: c.id)
def test_multiple_ov_blck_files(case, tmp_path):
    run_seed_search(case, tmp_path)


# ==========================================
# PYTEST STANDALONE TESTS 
# ==========================================

def test_inv(tmp_path):
    p = tmp_path / "a.yaml"
    p.write_text("""
reference: "tests/inputs/test01.fa"
max_tries: 100
variant_sets:
    - type: "INV"  
      number: 1
      length_ranges: [[3, 3]]
    """)

    sim = SVSimulator(config_path=str(p))
    sim.run()

    svs = sim.svs
    assert len(svs) == 1
    assert len(svs[0].operations) == 1
    assert (svs[0].operations[0].transform == Transform(transform_type=TransformType.INV, is_in_place=True, divergence_prob=0))
    assert (svs[0].operations[0].target_region.replace(order_key=()) == svs[0].operations[0].source_region)


def test_inv_exact(tmp_path):
    a_bed = tmp_path / "a.bed"
    a_bed.write_text("chr19\t0\t3\tLINE1\n")
    cfg = tmp_path / "a.yaml"
    cfg.write_text(f"""
reference: "tests/inputs/test01.fa"
max_tries: 1
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "INV"
      number: 1
      length_ranges: [[null, null]]
      overlap_region_type: [["LINE1"]]
      overlap_mode: exact
    """)

    sim = SVSimulator(config_path=str(cfg))
    sim.run()

    svs = sim.svs
    assert (svs[0].operations[0].target_region.replace(order_key=()) == Region(chrom='chr19', start=0, end=3))


def test_inv_contained(tmp_path):
    a_bed = tmp_path / "a.bed"
    a_bed.write_text("chr19\t3\t7\tLINE1\n")
    cfg = tmp_path / "a.yaml"
    cfg.write_text(f"""
reference: "tests/inputs/test01.fa"
max_tries: 1
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "INV"  
      number: 1
      length_ranges: [[3, 3]]
      overlap_region_type: [["LINE1"]]
      overlap_mode: contained
    """)

    sim = SVSimulator(config_path=str(cfg))
    sim.run()

    svs = sim.svs
    assert (svs[0].operations[0].target_region.replace(order_key=()) in
            [Region(chrom='chr19', start=3, end=6), Region(chrom='chr19', start=4, end=7)])


def test_inv_partial(tmp_path):
    a_bed = tmp_path / "a.bed"
    a_bed.write_text("chr19\t3\t7\tLINE1\n")
    cfg = tmp_path / "a.yaml"
    cfg.write_text(f"""
reference: "tests/inputs/test01.fa"
max_tries: 1
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "INV"  
      number: 1
      length_ranges: [[3, 3]]
      overlap_region_type: [["LINE1"]]
      overlap_mode: partial
    """)

    sim = SVSimulator(config_path=str(cfg))
    sim.run()

    svs = sim.svs
    assert (svs[0].operations[0].target_region.replace(order_key=()) in
            [Region(chrom='chr19', start=0, end=3),
             Region(chrom='chr19', start=1, end=4),
             Region(chrom='chr19', start=2, end=5),
             Region(chrom='chr19', start=5, end=8),
             Region(chrom='chr19', start=6, end=9)])


def test_trEXP(tmp_path):
    a_bed = tmp_path / "a.bed"
    a_bed.write_text("chrA\t4\t16\tALU\tTCG\n")
    cfg = tmp_path / "a.yaml"
    cfg.write_text(f"""
reference: "tests/inputs/test_tr.fa"
max_tries: 1
homozygous_only: true
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "trEXP"
      number: 1
      repeat_count_change_range: [2, 2]
      overlap_region_type: [["ALU"]]
    """)

    sim = SVSimulator(config_path=str(cfg))
    sim.run()

    sim_fa = tmp_path / "sim.hapA.fa"
    with FastaFile(str(sim_fa)) as fasta_file:
        hap = fasta_file.fetch(fasta_file.references[0])
        
    assert hap == 'AAAATCGTCGTCGTCGTCGTCGAAAA'


def test_trCON(tmp_path):
    a_bed = tmp_path / "a.bed"
    a_bed.write_text("chrA\t4\t16\tALU\tTCG\n")
    cfg = tmp_path / "a.yaml"
    cfg.write_text(f"""
reference: "tests/inputs/test_tr.fa"
max_tries: 1
homozygous_only: true
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "trCON"
      number: 1
      repeat_count_change_range: [3, 3]
      overlap_region_type: [["ALU"]]
    """)

    sim = SVSimulator(config_path=str(cfg))
    sim.run()

    sim_fa = tmp_path / "sim.hapA.fa"
    with FastaFile(str(sim_fa)) as fasta_file:
        hap = fasta_file.fetch(fasta_file.references[0])
        
    assert hap == 'AAAATCGAAAA'