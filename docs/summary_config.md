```yaml
# YAML CONFIG FILE
sim_settings:
    max_tries: 200
    prioritize_top: True
SVs:
    - type: "dDUP"  # "A_" -> "A_A"
      number: 3
      min_length:
        - 500  # minimum size for block A
        - 5000  # minimum size for block _
      max_length:
        - 1000  # minimum size for block A
        - 10000  # minimum size for block _
    - type: "INV_dDUP"  # "A_" -> "A_a'"
      number: 3
      min_length:
        - 500  # minimum size for block A
        - 5000  # minimum size for block _
      max_length:
        - 1000  # minimum size for block A
        - 10000  # minimum size for block _
    - type: "TRA_UNBALANCED"  # "A_" -> "_A"
      number: 3
      min_length:
        - 500  # minimum size for block A
        - 5000  # minimum size for block _
      max_length:
        - 1000  # maximum size for block A
        - 10000  # maximum size for block _
    - type: "TRA_BALANCED"  # "A_B" -> "B_A"
      number: 3
      min_length:
        - 500  # minimum size for block A
        - 5000  # minimum size for block _
        - 500  # minimum size for block B
      max_length:
        - 1000  # maximum size for block A
        - 10000  # maximum size for block _
        - 1000  # maximum size for block B
    - type: "DEL"  # "A" -> ""
      number: 3
      min_length: 1000  # minimum size for block A
      max_length: 10000  # maximum size for block A
    - type: "DUP"
      number: 3
      min_length: 1000
      max_length: 10000
    - type: "INV"
      number: 3
      min_length: 1000
      max_length: 10000
    - type: "INS"
      number: 3
      min_length: 1000
      max_length: 10000
    - type: "dupINVdup"  # "ABC" -> "Ac'ba'C"
      number: 3
      min_length:
        - 500  # minimum size for block A
        - 500  # minimum size for block B
        - 500  # minimum size for block C
      max_length:
        - 1000  # maximum size for block A
        - 1000  # maximum size for block B
        - 1000  # maximum size for block C
    - type: "delINVdup"
      number: 3
      min_length:
        - 500
        - 500
        - 500
      max_length:
        - 1000
        - 1000
        - 1000
    - type: "dupINVdel"
      number: 3
      min_length:
        - 500
        - 500
        - 500
      max_length:
        - 1000
        - 1000
        - 1000
    - type: "delINVdel"
      number: 3
      min_length:
        - 500
        - 500
        - 500
      max_length:
        - 1000
        - 1000
        - 1000
    - type: "INVdel"  # "AB" -> "b"
      number: 3
      min_length:
        - 500  # minimum size for block A
        - 500  # minimum size for block B
      max_length:
        - 1000  # maximum size for block A
        - 1000  # maximum size for block B
    - type: "delINV"
      number: 3
      min_length:
        - 500
        - 500
      max_length:
        - 1000
        - 1000
    - type: "dDUP_iDEL"  # "A_B" -> "A_A'"
      number: 3
      min_length:
        - 500  # minimum size for block A
        - 500  # minimum size for block B
        - 5000  # minimum size for block _
      max_length:
        - 1000  # maximum size for block A
        - 1000  # maximum size for block B
        - 10000  # maximum size for block _
    - type: "INS_iDEL"  # "A_B" -> "_A'"
      number: 3
      min_length:
        - 500
        - 500
        - 5000
      max_length:
        - 1000
        - 1000
        - 10000
    - type: "INVdup"  # "A" -> "aa'"
      number: 3
      min_length:
        - 500  # minimum size for block A
      max_length:
        - 1000  # maximum size for block A
    - type: "dup_INV"  # "AB" -> "Aba'"
      number: 3
      min_length:
        - 500  # minimum size for block A
        - 500  # minimum size for block B
      max_length:
        - 1000  # maximum size for block A
        - 1000  # maximum size for block B
    - type: "INV_dup"  # "AB" -> "b'aB"
      number: 3
      min_length:
        - 500
        - 500
      max_length:
        - 1000
        - 1000
    - type: "SNP"  # "A" -> "A*" (for A of length 1)
      number: 3
    - type: "DIVERGENCE"  # "A" -> "A*" (for A of arbitrary length)
      number: 3
      divergence_prob: 0.2  # optional parameter giving the probability that each base in the interval will be changed
      min_length:
        - 500
      max_length:
        - 1000
```
