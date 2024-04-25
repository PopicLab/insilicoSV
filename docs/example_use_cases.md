# Example Use Cases
insilicoSV provides various simulation features that can be used together or separately to generate synthetic genomes with varying levels of control over SV placement. Examples of the different use cases are provided below.

### Example 1 - Predefined SV types
To incorporate SVs from the predefined library of SV types, a configuration file of the following form can be provided with parameters provided for the count and min/max length for each set of SVs to be included in the output genome.
```yaml
# YAML config file
sim_settings:
    reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "INS"  # "" -> "A"
      number: 10
      length_ranges: [[5, 10]]
    - type: "INVdel"  # "AB" -> "a"
      number: 2
      length_ranges:
        - [5, 10]  # min/max length for INV
        - [10, 15]  # min/max length for DEL
    - type: "dupINVdel"  # "ABC" -> "Aba"
      number: 1
      length_ranges:
        - [5, 10]  # min/max length for first DEL
        - [10, 15]  # min/max length for INV
        - [5, 10]  # min/max length for second DEL
```

### *Example 1a* - Example config with entire insilicoSV vocabulary
This [summary config](summary_config.md) contains an example specifying the event size ranges for each predefined variant type,
as well as examples of the various placement features that can be used to bias variants towards or away from regions of interest.

### Example 1b - Example SNP/INDEL specification
SNPs and INDELs are supported by insilicoSV as well. Because SNPs are only a single base in length, they only need to be specified with `number` as in the example below:
```yaml
variant_sets:
  - type: "SNP"  # "A" -> "A*" (for A of length 1)
    number: 10
```
INDELs are not given a unique type label but can be simulated by setting a sufficiently small min and max size for an SV of type DEL or INS.

### Example 1c - Unbounded dispersions
When specifying length ranges for dispersion intervals, `None` can be provided as an upper bound and this will result in
the dispersion target locus being placed at any position on the same chromosome as the source event. For example:
```yaml
variant_sets:
  - type: "dDUP"  # "A_" -> "A_A"
    number: 2
    length_ranges:
      - [500, 1000]  # min/max length for source
      - [500, None]  # unbounded dispersion: target locus may be placed arbitrarily far
```
Unbounded dispersions can be used with predefined or custom SVs.

### Example 2 - Custom SVs
Custom SVs can be specified by manually describing the desired variant with the grammatical notation described in [SV grammar](sv_grammar.md). An example input config and output BED file are given below:
```yaml
# YAML config file
sim_settings:
    reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "Custom"
      source: "AB_C_D"
      target: "bb_AEc_EDC"
      number: 1
      length_ranges: 
        - [5, 10]   # A
        - [6, 10]   # B
        - [7, 10]   # C
        - [8, 10]   # D
        - [10, 15]  # E
        - [10, 15]  # first _
        - [15, 20]  # second _
```
```
# BED file
chr21	100	110	chr21	100	110	INV	    9	  0/1	AB_C_D>bb'_AEc'_EDC	1	0
chr21	100	110	chr21	109	110	INVDUP	9	  0/1	AB_C_D>bb'_AEc'_EDC	1	1
chr21	92	101	chr21	124	125	TRA	    8	  0/1	AB_C_D>bb'_AEc'_EDC	1	1   # order important for insertion-like operations at the same position
chr21	124	125	chr21	124	125	INS	    14	0/1	AB_C_D>bb'_AEc'_EDC	1	2
chr21	124	135	chr21	124	125	INVDUP	10	0/1	AB_C_D>bb'_AEc'_EDC	1	3
chr21	150	151	chr21	150	151	INS	    14	0/1	AB_C_D>bb'_AEc'_EDC	1	1
chr21	124	135	chr21	159	160	TRA	    10	0/1	AB_C_D>bb'_AEc'_EDC	1	1
```

### Example 3 - Editing reference with known SVs
To edit an input reference file with a known set of variants the user can provide them in a VCF file. 
The predefined variant types supported for this use case are DEL, DUP, INV, INS, dDUP, INV_dDUP, TRA_UNBALANCED, INVdup, 
and SNP. For insertions, VCF records may be specified with the insertion sequence given in an INFO field called `INSSEQ`
(provided a matching header line is included as well). All VCF records are expected to include an info field `SVTYPE`
to record variant type. The commandline call to perform this reference edit is:
```yaml
# YAML config file
sim_settings:
    reference: "{path}/{to}/ref.fa"
variant_sets:
    - vcf_path: "{path_to_vcf}"
```
```
insilicosv <config.yaml>
```

### Example 4 - Marking genome blacklist intervals
When initializing a new simulation the user can include a blacklist of genome intervals (i.e., intervals that specified
variants will avoid) via a VCF or BED file given (or multiple provided in a comma-separated list) in the `blacklist_regions`
entry of the config file:
```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
blacklist_regions: "{path}/{to}/{banned_intervals}.{bed, vcf}"
variant_sets:
    - type: "DEL"  # "A" -> ""
      number: 3
      length_ranges: [[1000, 10000]]
      blacklist_region_type: "all"
    ...
```
Each variant set specified in the config can include a `blacklist_region_type` value which will control which intervals
recorded in `blacklist_regions` will be avoided by those variants. In the above example, `"all"` is specified, which will
result in none of the three deletions from being placed in any of the regions in `blacklist_regions`. If `blacklist_regions`
is given input in BED file format, records can include a fourth column recording `type`, which can then be used to filter
the blacklist intervals considered for a given set of SVs.

If entries are provided in VCF format an arbitrary record ID and SVTYPE can be provided (e.g., 'EMPTY')

### Example 4b - Specifying a minimum inter-variant distance
Variant placement can also be constrained by enforcing that there be a minimum inter-variant distance between any two
breakpoints belonging to different variants. This value will be by default 1 (i.e., enforcing that no two variants
are separated by at least 1 bp) but a desired simulation-specific distance can be provided in the config input `min_intersv_dist`:
```yaml
sim_settings:
  reference: "{path}/{to}/ref.fa"
  min_intersv_dist: 100  # <- each DEL will be separated by at least 100bp
variant_sets:
  - type: "DEL"
    number: 3
    length_ranges: [(1000, 10000)]
```


### Example 5 - Placing SVs at known regions of interest
To augment a randomized simulation of SVs onto an input reference, the user can include in the simulation config
file the path to a BED file (or multiple) containing known element intervals (e.g., known repetitive elements taken from
RepeatMasker). For each variant set, if the `overlap_type` field is given a value of `"exact"` or `"partial"` then each
variant of that set will be placed in exact or partial overlap with a randomly selected interval from `overlap_regions`. 
An example config with these inputs is:
```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
overlap_regions: "/{path_to}/{candidate_overlap_events}.bed"
variant_sets:
    - type: "DEL"  # "A" -> ""
      number: 5
      length_ranges: [[5, 5]]
      overlap_type: "exact"  # <- or "partial"
```
Multiple BED files can be given as input and their records will be combined and drawn from during SV placement (in this
case the user should provide a list of paths). Events from the `overlap_regions` BED file(s) will be shuffled on input,
and the file is required to have the first four columns of standard BED records (chrom, chromStart, chromEnd, name).

For each variant set the user can specify overlap with different interval types using the `overlap_element` field. If
an `overlap_element` value is provided, the intervals extracted from `overlap_regions` will be limited to those with a
type (recorded in the fourth BED record column) with the `overlap_element` value as a prefix. For example, the below
config will yield two deletions placed completely at random, three placed in exact match with L1HS intervals, and
two placed in partial overlap with L1PA3 intervals:
```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
overlap_regions: ['/{path_to}/{candidate_overlap_events_1}.bed','/{path_to}/{candidate_overlap_events_2}.bed']
variant_sets:
    - type: "DEL"  # "A" -> ""
      number: 2
      length_ranges: [[500, 1000]]
    - type: "DEL"  # "A" -> ""
      number: 3
      length_ranges: [[500, 1000]]
      overlap_type: "exact"
      overlap_element: "L1HS"
    - type: "DEL"  # "A" -> ""
      number: 2
      length_ranges: [[500, 1000]]
      overlap_type: "partial"
      overlap_element: "L1PA3"
```
The output VCF file will label which SVs were placed at specified intervals with the additional INFO field
`OVERLAP_EV={evt. name}', as in this example record:
```
chr21   18870078    DEL N   DEL 100 PASS    END=18876908;SVTYPE=DEL;SVLEN=6831;OVERLAP_EV=L1HS  GT  0/1
```

### Example 5a - Placing specific SV components at regions of interest
In the case of complex SVs that include multiple sub-SV events, there are various ways in which the user can specify 
which SV component will overlap a region of interest. For SVs that are contiguous (i.e., don't include a dispersion)
overlap can be assigned to a random or specific event within the SV, or it can be assigned to the full SV interval.
Random or full-SV overlap can be specified using the `overlap_component` field as in the example below:
```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
overlap_regions: "/{path_to}/{candidate_overlap_events_1}.bed"
variant_sets:
    - type: "delINVdel"  # "ABC" -> "b"
      number: 5
      length_ranges:
        - [500, 1000]
        - [500, 1000]
        - [500, 1000]
      overlap_type: "exact"
      overlap_region_type: "L1HS"
      overlap_component: "rand"  # <- or "full_sv"
```
To specify overlap with a specific event, `overlap_region_type` can take a list of length equal to the number of events
there are in the SV. In the below example, each delINVdel will be placed such that the middle INV event exactly overlaps
an L1HS drawn from `overlap_regions`.
```yaml
sim_settings:
  reference: "{path}/{to}/ref.fa"
overlap_regions: "/{path_to}/{candidate_overlap_events}.bed"
variant_sets:
  - type: "delINVdel"  # "ABC" -> "b"
    number: 5
    length_ranges:
      - [500, 1000]
      - [500, 1000]
      - [500, 1000]
    overlap_type: "exact"
    overlap_region_type: [None, "L1HS", None]
```
In the case of providing `overlap_region_type` with a list, it is required that only one entry be non-None in order to
ensure that the SV can be placed properly. This feature can be used in the same way with custom SVs, as in the example below:
```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
overlap_events: '/{path_to}/{candidate_overlap_events}.bed'
variant_sets:
    - type: "Custom"
      source: "AB_C_D"
      target: "bb'_AEc'_EDC"
      number: 1
      length: 
        - [5, 10]   # A
        - [6, 10]   # B
        - [7, 10]   # C
        - [8, 10]   # D
        - [10, 15]  # E
        - [10, 15]  # first _
        - [15, 20]  # second _
      overlap_type: "partial"
      overlap_region_type: [None, "L1HS", None, None, None, None, None]
```

For SVs involving dispersions, overlap can be assigned to the source or target locus also using the `overlap_component`
field by providing values of `source` or `target`. E.g.:
```yaml
sim_settings:
  reference: "{path}/{to}/ref.fa"
overlap_regions: "/{path_to}/{candidate_overlap_events_1}.bed"
variant_sets:
  - type: "dDUP"  # "A_" -> "A_A"
    number: 5
    length_ranges:
      - [500, 1000]
      - [500, 1000]
    overlap_type: "exact"
    overlap_component: "source"  # <- or "target"
```
If target overlap is specified, the SV's target locus is placed at a random point within the selected interval. In this
case `overlap_type` must be left blank (as there is no notion of partial or exact overlap for the target locus).

Additionally, `overlap_region_type` can also be provided a list for specific event overlap as in the previous example,
with the additional option that multiple events can be assigned overlaps if they are separated by an unbounded dispersion.
For example, the below config will yield a dDUP with both source and target being placed in L1HS intervals, but this
is only possible because the dispersion is unbounded and can therefore be set to whatever size is necessary to admit
L1HS placement for both events:
```yaml
sim_settings:
  reference: "{path}/{to}/ref.fa"
overlap_regions: "/{path_to}/{candidate_overlap_events}.bed"
variant_sets:
  - type: "dDUP"  # "A_" -> "A_A"
    number: 5
    length_ranges:
      - [500, 1000]
      - [500, None]
    overlap_type: "exact"
    overlap_region_type: ["L1HS", "L1HS"]
```


### Example 5a - Placing DUPs or DELs at Alu-mediated intervals
An additional use case for the above known-element placement is to place deletion or tandem duplication SVs in between
Alu elements in the genome (Alu-mediated DEL/DUPs being a well-studied case of SV/repetitive element relation – e.g., 
[Gu et al., 2015](https://academic.oup.com/hmg/article/24/14/4061/2385874)). Alu-mediated DELs or DUPs can be specified 
in the same way as the above cases of specifying overlap, but instead with the inputs `overlap_type: "flanked"` and
`overlap_region_type: "Alu"`. An example config file is given below:
```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
overlap_regions: "/{path_to}/{candidate_overlap_events}.bed"
variant_sets:
    - type: "DEL"  # "A" -> ""
      number: 10
      length_ranges: [[500, 1000]]
      overlap_type: "flanked"
      overlap_region_type: "Alu"
```

### Example 6 - Divergent intervals
In addition to the various SV types included in the predefined library, insilicoSV also allows for the simulation of divergent
intervals in which a random proportion of the bases in a given interval will be corrupted with some probability `p`. Divergent
intervals can be included in the same way as any of the SVs provided, accessible through the built-in type name `'DIVERGENCE'` (or through custom SVs using the asterisk divergence operator as described in [SV grammar](sv_grammar.md)). The
probability parameter by which each base will be randomly changed can be optionally provided by the user as shown in the example
below (and if it is not provided it will be drawn uniformly from (0.5, 1.0)):
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
variant_sets:
    - type: "DIVERGENCE"  # "A" -> "A*"
      number: 3
      divergence_prob: 0.2
      length_ranges: [[500, 1000]]
```
