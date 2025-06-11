```yaml
# YAML CONFIG FILE
reference: "{path}/{to}/ref.fa"
overlap_regions: ["{path}/{to}/{candidate_overlap_events_1}.bed"]
blacklist_regions: ["{path}/{to}/{candidate_overlap_events_2}.bed"]
variant_sets:
    # ==== PREDEFINED VARIANT EXAMPLES ====
    - type: "DEL"  # "A" -> ""
      number: 3
      length_ranges: [[1000, 10000]]  # minimum, maximum size for A
    - type: "DUP"  # "A" -> "AA+" by default n_copies is [1] for DUP
      number: 3
      length_ranges: [[1000, 10000]]
    - type: "mCNV"  # "A" -> "A+"
      number: 3
      n_copies: [5]
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
    - type: "INV_dDUP"  # "A_" -> "A_a"
      number: 3
      length_ranges:
        - [500, 1000]
        - [5000, 10000]
    - type: "nrTRA"  # "A_" -> "_A"
      number: 3
      length_ranges:
        - [500, 1000]
        - [5000, 10000]
    - type: "rTRA"  # "A_B" -> "B_A"
      number: 3
      length_ranges:
        - [500, 1000]  # minimum, maximum size for A
        - [500, 1000]  # minimum, maximum size for B
        - [5000, 10000]    # minimum, maximum size for _
    - type: "INVdel"  # "AB" -> "a"
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
    - type: "dDUP_iDEL"  # "A_B" -> "A_A"
      number: 3
      length_ranges:
        - [500, 1000] # minimum, maximum size for A
        - [500, 1000] # minimum, maximum size for B
        - [5000, 10000] # minimum, maximum size for _
    - type: "INS_iDEL"  # "A_B" -> "_A"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
        - [5000, 10000]
    - type: "DUP_INV"  # "A" -> "aa"
      number: 3
      length_ranges: [[500, 1000]]
    - type: "dupINV"  # "AB" -> "Aba"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
    - type: "INVdup"  # "AB" -> "baB"
      number: 3
      length_ranges:
        - [500, 1000]
        - [500, 1000]
    - type: "SNP"  
      number: 3
    # If null given as dispersion max length, target can be placed anywhere on the chromosome 
    # null is only allowed for dispersions and if the min length is null then the max_length has to null as well. 
    - type: "dDUP"  # "A_" -> "A_A"
      number: 5
      length_ranges:
        - [50, 100]
        - [500, null]
        
  # ==== OVERLAP PLACEMENT EXAMPLES ====
  # 1) 'overlap_mode' will result in each of these DELs being placed to overlap with a randomly
  #    chosen interval from 'overlap_regions'.  If 'overlap_mode' is 'exact',
  #    the DEL will exactly match the ROI.  If 'overlap_mode' is 'partial', the DEL
  #    will partially (but not fully) overlap the ROI.  If 'overlap_mode' is 'contained',
  #    the DEL will be fully contained within the ROI.
  #
  #    The ROIs to be chosen from can be limited by filtering conditions.
  #    "overlap_region_type" specifies region types to match.  Region type is given in column 4 of
  #    the .bed file.  A string from the "overlap_region_type" list will match .bed regions whose
  #    type begins with that string; the string "ALL" will match all region types.
  #
    - type: "DEL"  # "A" -> ""
      number: 5
      length_ranges: [[50, 100]]
      overlap_mode: "partial"  # or "contained" or "exact"
      overlap_region_type: ["L1HS"]  # optional filter on the type of intervals selected for overlap
  # 2) The SV part involved in the overlap is termed the overlap anchor.  For SVs without dispersions,
  # the specification of the overlap anchor can be omitted, in which case the anchor will default
  # to the full SV.
    - type: "delINVdel"  # "(ABC)" -> "b"
      number: 5
      length_ranges:
        - [50, 100]
        - [50, 100]
        - [50, 100]
      overlap_mode: "partial"
      overlap_region_type: ["L1"]
  # 3) The overlap anchor can be specified on either the source or target of the SV's grammar
  # definition, by putting in parentheses the part of the SV constituting the anchor.
    - type: "(A)_ -> A_A"   # Also allowed: "A_()" "(A_)" 
      number: 5
      length_ranges:
        - [50, 100]
        - [500, 1000]
      overlap_mode: "contained"


  # ==== BLACKLIST EXAMPLE ====
    - type: "DEL"  # "(A)" -> ""
      number: 5
      length_ranges: [[50, 100]]
      blacklist_region_type: ["all"]  # DEL placement will avoid all intervals listed in `blacklist_regions`
                                    # -> can alternatively filter by region type as in the above examples
```
