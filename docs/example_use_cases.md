# Example Use Cases
insilicoSV provides various simulation features that can be used together or separately to generate synthetic genomes with varying levels of control over SV placement. Examples of the different use cases are provided below.

### Example 1 - Predefined SVs
To incorporate SVs from the predefined library of SV types, a configuration file of the following form can be provided with parameters provided for the count and min/max length for each set of SVs to be included in the output genome.
```yaml
# YAML config file
sim_settings:
    reference: {path}/{to}/ref.fa
variant_sets:
    - type: "INS"
      number: 10
      min_length: [5]
      max_length: [10]
    - type: "INVdel"
      number: 2
      min_length: [5]
      max_length: [10]
    - type: "dupINVdel"
      number: 1
      min_length:
        - 5
        - 10
        - 5
      max_length:
        - 10
        - 15
        - 10
```

### *Example 1a* - Example config with entire insilicoSV vocabulary
This [summary config](summary_config.md) contains an example specifying the event size ranges for each predefined SV.

### Example 1b - Example SNP/INDEL specification
SNPs and INDELs are supported by insilicoSV as well. Because SNPs are only a single base in length, they only need to be specified with `number` as in the example below:
```yaml
variant_sets:
  - type: "SNP"
    number: 10
```
INDELs are not given a unique type label but can be simulated by setting a sufficiently small min and max size for an SV of type DEL or INS.

### Example 2 - Custom SVs
Custom SVs can be specified by manually describing the desired variant with the grammatical notation described in [SV grammar](sv_grammar.md). An example input config and output BED file are given below:
```yaml
# YAML config file
sim_settings:
    reference: {path}/{to}/ref.fa
variant_sets:
    - type: "Custom"
      source: AB_C_D
      target: bb'_AEc'_EDC
      number: 1
      min_length: 
        - 5   # A
        - 6   # B
        - 7   # C
        - 8   # D
        - 10  # E
        - 10  # first _
        - 15  # second _
      max_length: 
        - 10
        - 10
        - 10
        - 10
        - 15
        - 15
        - 20
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
To edit an input reference file with a known set of variants the user can provide them in a VCF file. The events in the VCF must be non-overlapping. All single-interval and dispersion-based predefined variant types are supported for this use case (i.e., DEL, DUP, INV, INS, dDUP, INV_dDUP, TRA, INVdup, and SNP). For insertions, VCF records may be specified with the insertion sequence given in an INFO field called `INSSEQ` (provided a matching header line is included as well). All VCF records are expected to include an info field `SVTYPE` to record variant type. The commandline call to perform this reference edit is:
```yaml
# YAML config file
sim_settings:
    reference: {path}/{to}/ref.fa
variant_sets:
    - vcf_path: {path_to_vcf}
```
```
insilicosv <config.yaml>
```

### Example 4 - Marking banned intervals of the genome
When initializing a new simulation the user can include a list of banned genome intervals (i.e., intervals in which no variants will be simulated) via a VCF given in the `avoid_intervals` entry of the `SVs` section of the config file:
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
avoid_intervals: "{path}/{to}/{banned_intervals}.vcf"
variant_sets:
    - type: "DEL"
      number: 3
      min_length: [1000]
      max_length: [10000]
    - type: "DUP"
      number: 3
    ...
```
The entries of the VCF will only have the interval information extracted under this feature, so an arbitrary record ID and SVTYPE can be provided (e.g., 'EMPTY')


### Example 5 - Placing SVs at known repetitive element intervals
To augment a randomized simulation of SVs onto an input reference, the user can include in the simulation config
file the path to a BED file (or multiple) containing known element intervals (e.g., known repetitive elements taken from RepeatMasker).
In addition to providing the path to the relevant BED file(s), the user will also need to specify how many of each 
SV type they wish to be placed at events from the BED file. An example config with these inputs is:
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
overlap_events:
    bed: '/{path_to}/{candidate_overlap_events}.bed'
variant_sets:
    - type: "DEL"
      number: 10
      min_length: [5]
      max_length: [5]
      num_overlap: 5
    - type: "DUP"
      number: 10
      min_length: [5]
      max_length: [5]
      num_overlap: 2
```
Multiple BED files can be given as input and their records will be combined and drawn from during SV placement (in this
case the user should provide a list of paths). Additionally, the user can provide a list of repetitive element types that
will be allowed for SV placement during simulation (BED records of all other types being ignored). An example entry
with multiple BED files and specified allowed types is:
```yaml
overlap_events:
    bed: ['/{path_to}/{candidate_overlap_events_1}.bed','/{path_to}/{candidate_overlap_events_2}.bed']
    allow_types: ['L1HS', 'L1PA3']
```
While the simulator is placing each set of SVs, the first (`num_overlap`) SVs of that type will have their location
given by an event from the BED file given in the `overlap_events` field that falls within the specified size range for 
that SV type. Events from the `overlap_events` BED file(s) will be shuffled on input, and the file is required to have
the first four columns of standard BED records (chrom, chromStart, chromEnd, name).

The labels in the `allow_types` field can be given either as full element names or prefixes. For instance, if 'L1' is provided
then all elements in the input BED file(s) with a name beginning with 'L1' will be considered for selection.

The output VCF file will label which SVs were placed at specified intervals with the additional INFO field
`OVERLAP_EV={evt. name}', as in this example record:
```
chr21   18870078    DEL N   DEL 100 PASS    END=18876908;SVTYPE=DEL;SVLEN=6831;OVERLAP_EV=L1HS  GT  0/1
```
For SVs involving dispersions (dDUP, INV_dDUP, TRA) the position is assigned such that the source event of the SV (the component
getting duplicated or translocated) is placed at the selected element interval. For complex SVs with multiple non-dispersion
source fragments (e.g., delINVdel), one of the non-dispersion source fragments is chosen at random to be the overlapping
component with the selected known element.

### Example 5a - Placing DUPs or DELs at Alu-mediated intervals
An additional use case for the above known-element placement is to place deletion or tandem duplication SVs in between Alu elements in the genome (Alu-mediated CNVs being a well-studied case of SV/repetitive element relation – e.g., [Gu et al., 2015](https://academic.oup.com/hmg/article/24/14/4061/2385874)). Alu-mediated DELs or DUPs can be specified in the same way as the above cases of specifying overlap, but instead by specifying the desired number of Alu-mediated SVs with the config field `num_alu_mediated`. An example config file is given below:
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
overlap_events:
    bed: '/{path_to}/{candidate_overlap_events}.bed'
variant_sets:
    - type: "DEL"
      number: 10
      min_length: [500]
      max_length: [1000]
      num_alu_mediated: 5
```

### Example 5b - Specifying different overlap counts for different element types
In situations where a list of different `allow_types` are given alongside the input BED file of known elements, different counts can be given for the different element types listed in the `allow_types` field. For example, the augmented version of the above example shows how one might specify that a different number of DELs should be places at L1HS intervals than should be placed at L1PA3 intervals:
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
overlap_events:
    bed: ['/{path_to}/{candidate_overlap_events_1}.bed','/{path_to}/{candidate_overlap_events_2}.bed']
    allow_types: ['L1HS', 'L1PA3']
variant_sets:
    - type: "DEL"
      number: 10
      min_length: [500]
      max_length: [1000]
      num_overlap: [3, 5]
```
By providing a list in the `num_overlap` field, each number given in the list will be interpreted as the desired overlap count for the corresponding entry in the `allow_types` list. As a result, any list given in the `num_overlap` field must be of the same length as the `allow_types` list. `0` is a valid entry in the `num_overlap` list, and are required if one wishes to only specify overlap counts for a subset of the element types given in `allow_types`.

### Example 5c - Partial Overlap
In the same way that SVs can be made to completely overlap known element intervals in the various ways described above, they can also be made to *partially* overlap known element intervals. This can be done using SV config field `num_partial_overlap` which operates in the same way as `num_overlap`. An analogous version of the above example with partial overlap is given below:
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
overlap_events:
    bed: ['/{path_to}/{candidate_overlap_events_1}.bed','/{path_to}/{candidate_overlap_events_2}.bed']
    allow_types: ['L1HS', 'L1PA3']
variant_sets:
    - type: "DEL"
      number: 10
      min_length: [500]
      max_length: [1000]
      num_partial_overlap: [3, 5]
```

**The various types of known element placement are not mutually exclusive and can be used concurrently (when appropriate with respect to the types of SVs being simulated)**

### Example 6 - Divergent intervals
In addition to the various SV types included in the predefined library, insilicoSV also allows for the simulation of divergent
intervals in which a random proportion of the bases in a given interval will be corrupted with some probability `p`. Divergent
intervals can be included in the same way as any of the SVs provided, accessible through the type name `'DIVERGENCE'`. The
probability parameter by which each base will be randomly changed can be optionally provided by the user as shown in the example
below (and if it is not provided it will be drawn uniformly from (0.5, 1.0)):
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
variant_sets:
    - type: "DIVERGENCE"
      number: 3
      divergence_prob: 0.2
      min_length: [500]
      max_length: [1000]
```
