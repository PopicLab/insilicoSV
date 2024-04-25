```yaml
# YAML CONFIG FILE
sim_settings:
    reference: "{path}/{to}/ref.fa"
overlap_regions: "{path}/{to}/{candidate_overlap_events_1}.bed"
blacklist_regions: "{path}/{to}/{candidate_overlap_events_2}.bed"
variant_sets:
    # ==== PREDEFINED VARIANT EXAMPLES ====
    - type: "DEL"  # "A" -> ""
      number: 3
      length_ranges: [[1000, 10000]]  # minimum, maximum size for A
    - type: "DUP"  # "A" -> "AA"
      number: 3
      length_ranges: [[1000, 10000]]
    - type: "INV"  # "A" -> "a"
      number: 3
      length_ranges: [[1000, 10000]]
    - type: "INS"  # "" -> "A"
      number: 3
      length_ranges: [[1000, 10000]]
    - type: "dDUP"  # "A_" -> "A_A"
      number: 3
      length_ranges:
        - [500, 1000]  # minimum, maximum size for A
        - [5000, 10000]  # minimum, maximum size for _
    - type: "INV_dDUP"  # "A_" -> "A_a'"
      number: 3
      length_ranges:
        - [500, 1000]
        - [5000, 10000]
    - type: "TRA_UNBALANCED"  # "A_" -> "_A"
      number: 3
      length_ranges:
        - [500, 1000]
        - [5000, 10000]
    - type: "TRA_BALANCED"  # "A_B" -> "B_A"
      number: 3
      length_ranges:
        - [500, 1000]  # minimum, maximum size for A
        - [500, 1000]  # minimum, maximum size for B
        - [5000, 10000]    # minimum, maximum size for _
    - type: "INVdel"  # "AB" -> "b"
      number: 3
      length_ranges:
        - [500, 1000] # minimum, maximum size for A
        - [500, 1000] # minimum, maximum size for B
    - type: "delINV"  # "AB" -> "b"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
    - type: "dupINVdup"  # "ABC" -> "AcbaC"
      number: 3
      length_ranges:
        - [500, 1000] # minimum, maximum size for A
        - [500, 1000] # minimum, maximum size for B
        - [500, 1000] # minimum, maximum size for C
    - type: "delINVdup"  # "ABC" -> "cbC"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
        - [500, 1000]
    - type: "dupINVdel"  # "ABC" -> "Aba"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
        - [500, 1000]
    - type: "delINVdel"  # "ABC" -> "b"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
        - [500, 1000]
    - type: "dDUP_iDEL"  # "A_B" -> "A_A'"
      number: 3
      length_ranges:
        - [500, 1000] # minimum, maximum size for A
        - [500, 1000] # minimum, maximum size for B
        - [5000, 10000] # minimum, maximum size for _
    - type: "INS_iDEL"  # "A_B" -> "_A'"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
        - [5000, 10000]
    - type: "INVdup"  # "A" -> "aa'"
      number: 3
      length_ranges: [[500, 1000]]
    - type: "dup_INV"  # "AB" -> "Aba'"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
    - type: "INV_dup"  # "AB" -> "b'aB"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
    - type: "SNP"  # "A" -> "A*" (for A of length 1)
      number: 3
    - type: "DIVERGENCE"  # "A" -> "A*" (for A of arbitrary length)
      number: 3
      divergence_prob: 0.2  # optional parameter giving the probability that each base in the interval will be changed
      length_ranges: [[500, 1000]]
  # ==== OVERLAP PLACEMENT EXAMPLES ====
  # 1) 'overlap_type' will result in each of these DELs being placed in exact (or partial) overlap with a randomly
  #    chosen interval from 'overlap_regions' of appropriate size
    - type: "DEL"  # "A" -> ""
    - number: 5
    - length_ranges: [[50, 100]]
    - overlap_type: "exact"  # or "partial"
    - overlap_region_type: "L1HS"  # optional filter on the type of intervals selected for overlap
  # 2a) If 'overlap_region_type' is given as a list, it will be interpreted as the fragment-level assignment of overlap
    - type: "delINVdel"  # "ABC" -> "b"
    - number: 5
    - length_ranges: 
        - [50, 100]
        - [50, 100]
        - [50, 100]
    - overlap_type: "exact"
    - overlap_region_type: ["L1HS", None, None]  # will result in the first DEL overlapping exactly with an L1HS
  # 2b) For SVs not involving a dispersion, can choose a random fragment or the full SV interval for overlap 
    - type: "delINVdel"  # "ABC" -> "b"
    - number: 5
    - length_ranges:
        - [50, 100]
        - [50, 100]
        - [50, 100]
    - overlap_type: "exact"
    - overlap_region_type: "L1HS"
    - overlap_component: "rand"  # or "full_sv"
  # 3) For SVs involving a dispersion, can assign overlap to the source or target interval with "source" or "target" 
    - type: "dDUP"  # "A_" -> "A_A"
    - number: 5
    - length_ranges:
        - [50, 100]
        - [500, 1000]
    - overlap_type: "exact"
    - overlap_component: "source"  # or "target" (for "target" overlap placing the target locus at a random point
    #                                in the selected overlap interval; if "target" given, must omit overlap_type)
  # 4) If None given as dispersion max length, target can be placed anywhere on the chromosome 
    - type: "dDUP"  # "A_" -> "A_A"
    - number: 5
    - length_ranges:
        - [50, 100]
        - [500, None]
    - overlap_type: "exact"
    # - overlap_component: ["L1HS", "SVA"]  <- and if dispersion is unbounded, can place both source and target
    #                                          at overlap positions
  # ==== BLACKLIST EXAMPLE ====
    - type: "DEL"  # "A" -> ""
    - number: 5
    - length_ranges: [[50, 100]]
    - blacklist_region_type: "all"  # DEL placement will avoid all intervals listed in `blacklist_regions`
                                    # -> can alternatively filter by region type as in the above examples
```
