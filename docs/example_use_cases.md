# Example Use Cases
insilicoSV provides various simulation modes that can be used together or separately to generate synthetic genomes with varying levels of control over SV placement. Examples of the different use cases are provided below.

### Example 1 - Predefined SVs
To incorporate SVs from the predefined library of SV types, a configuration file of the following form can be provided with parameters provided for the count and min/max length for each set of SVs to be included in the output genome.
```yaml
# YAML config file
sim_settings:
    max_tries: 200
    prioritize_top: True
SVs:
    - type: "INS"
      number: 10
      min_length: 5
      max_length: 10
    - type: "INVdel"
      number: 2
      min_length: 5
      max_length: 10
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
The output of the simulation run with the above config file would include a `.bed` file such as the following describing the placement of each SV:
```
# BED file
Chromosome21	148	149 Chromosome21	148	149	INS	10	  0/1	INS	        1	  1
Chromosome19	5	  6	  Chromosome19	5	  6	  INS	6	    0/1	INS	        2	  1
Chromosome19	38	39	Chromosome19	38	39	INS	7	    1/1	INS	        3   1
Chromosome21	48	49	Chromosome21	48	49	INS	8	    1/1	INS	        4	  1
Chromosome19	86	87	Chromosome19	86	87	INS	10	  0/1	INS	        5	  1
Chromosome19	64	65	Chromosome19	64	65	INS	9	    1/1	INS	        6	  1
Chromosome19	7	  8	  Chromosome19	7	  8	  INS	10	  0/1	INS	        7	  1
Chromosome21	141	142	Chromosome21	141	142	INS	8	    1/1	INS	        8	  1
Chromosome19	74	75	Chromosome19	74	75	INS	10	  1/1	INS	        9	  1
Chromosome19	60	61	Chromosome19	60	61	INS	7	    0/1	INS	        10  1
Chromosome21	23	31	Chromosome21	23	31	INV	7	    0/1	INVdel	    11	0
Chromosome21	30	36	Chromosome21	30	31	DEL	5	    0/1	INVdel	    11	0
Chromosome21	122	132	Chromosome21	122	132	INV	9	    1/1	INVdel	    12	0
Chromosome21	131	142	Chromosome21	131	132	DEL	10	  1/1	INVdel	    12	0
Chromosome19	93	106	Chromosome19	93	106	INV	12	  0/1	dupINVdel	  13	0
Chromosome19	88	94	Chromosome19	105	106	INVDUP	5	0/1	dupINVdel	  13	1
Chromosome19	105	113	Chromosome19	105	106	DEL	7	    0/1	dupINVdel	  13	0
```
### *Example 1a* - Example config with entire insilicoSV vocabulary
This [gist](https://gist.github.com/crohlicek/9d529e600508870b1424d1f41215acb8) contains an example config file specifying the subevent size ranges for each event.

### Example 1b - Example SNP specification
We include SNPs as an available event type for simulation and admit them to the input config files with a modified form of the default SV config info used for other events.
Because SNPs are only a single base in length, they only need to be specified with `number` as in the example below:
```yaml
SVs:
  - type: "SNP"
    number: 10
```
Usage of this event type is otherwise the same.

### Example 2 - Custom SVs
Similar to specifying predefined SVs for a simulation, custom SVs can be specified by manually describing the desired variant with the grammatical notation described in [SV grammar](sv_grammar.md). An example input config and output `.bed` file are given below:
```yaml
# YAML config file
SVs:
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
# BEDPE file
Chromosome21	100	110	Chromosome21	100	110	INV	    9	  0/1	AB_C_D>bb'_AEc'_EDC	1	0
Chromosome21	100	110	Chromosome21	109	110	INVDUP	9	  0/1	AB_C_D>bb'_AEc'_EDC	1	1
Chromosome21	92	101	Chromosome21	124	125	TRA	    8	  0/1	AB_C_D>bb'_AEc'_EDC	1	1   # order important for insertion-like operations at the same position
Chromosome21	124	125	Chromosome21	124	125	INS	    14	0/1	AB_C_D>bb'_AEc'_EDC	1	2
Chromosome21	124	135	Chromosome21	124	125	INVDUP	10	0/1	AB_C_D>bb'_AEc'_EDC	1	3
Chromosome21	150	151	Chromosome21	150	151	INS	    14	0/1	AB_C_D>bb'_AEc'_EDC	1	1
Chromosome21	124	135	Chromosome21	159	160	TRA	    10	0/1	AB_C_D>bb'_AEc'_EDC	1	1
```

### Example 3 - Editing reference with input SVs
To edit an input reference file with a known set of SVs the user can provide a VCF file containing the SVs in the yaml 
of format shown above. The events in the VCF must be non-overlapping. All single-interval and dispersion-based predefined variant types are supported for this use case
(i.e., DEL, DUP, INV, INS, dDUP, INV_dDUP, TRA, INVdup, and SNP). For insertions, events may be specified with the insertion
sequence given in an INFO field called `INSSEQ` (provided a matching header line is included as well). All VCF records are
expected to include an info field `SVTYPE` to record event type. The commandline call to perform this reference edit is:
```yaml
# YAML config file
SVs:
    - vcf_path: {path_to_vcf}
```
```
insilicosv <ref.fna> <par.yaml> <output_prefix>
```

### Example 4 - Marking banned intervals of the genome
When initializing a new simulation the user can include a list of banned genome intervals (i.e., intervals in which no SVs will be simulated) via a VCF given in the `avoid_intervals` entry of the `SVs` section of the config file:
```yaml
sim_settings:
    max_tries: 200
    prioritize_top: True
SVs:
    - avoid_intervals: "{path}/{to}/{banned_intervals}.vcf"
    - type: "DEL"
      number: 3
      min_length: 1000
      max_length: 10000
    - type: "DUP"
      number: 3
    ...
```
The entries of the VCF will only have the interval information extracted under this feature, so an arbitrary record ID and SVTYPE can be provided (e.g., 'EMPTY')


### Example 5 - Placing events at known repetitive element intervals
To augment a randomized simulation of events onto an input reference, the user can include in the simulation config
file the path to a .bed file containing known element intervals (e.g., known repetitive elements taken from RepeatMasker).
In addition to providing the path to the relevant .bed file(s), the user will also need to specify how many of each 
event type they wish to be placed at events from the .bed file. An example config with these inputs is:
```yaml
sim_settings:
    max_tries: 200
    prioritize_top: True
overlap_events:
    bed: '/{path_to}/{candidate_overlap_events}.bed'
SVs:
    - type: "DEL"
      number: 10
      min_length:
        - 5
      max_length:
        - 5
      num_overlap: 5
    - type: "DUP"
      number: 10
      min_length:
        - 5
      max_length:
        - 5
      num_overlap: 2
```
Multiple .bed files can be given as input and their records will be combined and drawn from during event placement (in this
case the user should provide a list of paths). Additionally, the user can provide a list of repetitive element types that
will be allowed for event placement during simulation (.bed records of all other types being ignored). An example entry
with multiple .bed files and specified allowed types is:
```yaml
overlap_events:
    bed: ['/{path_to}/{candidate_overlap_events_1}.bed','/{path_to}/{candidate_overlap_events_2}.bed']
    allow_types: ['L1HS', 'L1PA3']
```
While the simulator is placing each set of events, the first (`num_overlap`) events of that type will have their location
given by an event from the .bed file given in the `overlap_events` field that falls within the specified size range for 
that SV type. Events from the `overlap_events` .bed file(s) will be shuffled on input, and the file is required to have
the first four columns of standard .bed records (chrom, chromStart, chromEnd, name).

The labels in the `allow_types` field can be given either as full element names or prefixes. For instance, if 'L1' is provided
then all elements in the input .bed file(s) with 'L1' in the name will be considered for selection.

The output .vcf file will label which events were placed at specified intervals with the additional INFO field
`OVERLAP_EV={evt. name}', as in this example record:
```
chr21   18870078    DEL N   DEL 100 PASS    END=18876908;SVTYPE=DEL;SVLEN=6831;OVERLAP_EV=L1HS  GT  0/1
```
For events involving dispersions (dDUP, INV_dDUP, TRA) the position is assigned such that the source event of the SV (the component
getting duplicated or translocated) is placed at the selected element interval. For complex events with multiple non-dispersion
source fragments (e.g., delINVdel), one of the non-dispersion source fragments is chosen at random to be the overlapping
component with the selected known element.

### Example 5a - Placing DUPs or DELs at Alu-mediated intervals
An additional use case for the above known-element placement is to place deletion or tandem duplication events in between Alu elements in the genome (Alu-mediated CNVs being a well-studied case of SV/repetitive element relation – e.g., [Gu et al., 2015](https://academic.oup.com/hmg/article/24/14/4061/2385874)). Alu-mediated DELs or DUPs can be specified in the same way as the above cases of specifying overlap, but instead by specifying the desired number of Alu-mediated SVs with the config field `num_alu_mediated`. An example config file is given below:
```yaml
sim_settings:
    max_tries: 200
    prioritize_top: True
overlap_events:
    bed: '/{path_to}/{candidate_overlap_events}.bed'
SVs:
    - type: "DEL"
      number: 10
      min_length: 500
      max_length: 1000
      num_alu_mediated: 5
```

### Example 5b - Specifying different overlap counts for different element types
In situations where a list of different `allow_types` are given alongside the input .bed file of known elements, different counts can be given for the different element types listed in the `allow_types` field. For example, the augmented version of the above example shows how one might specify that a different number of DELs should be places at L1HS intervals than should be placed at L1PA3 intervals:
```yaml
overlap_events:
    bed: ['/{path_to}/{candidate_overlap_events_1}.bed','/{path_to}/{candidate_overlap_events_2}.bed']
    allow_types: ['L1HS', 'L1PA3']
SVs:
    - type: "DEL"
      number: 10
      min_length: 500
      max_length: 1000
      num_overlap: [3, 5]
```
By providing a list in the `num_overlap` field, each number given in the list will be interpreted as the desired overlap count for the corresponding entry in the `allow_types` list. As a result, any list given in the `num_overlap` field must be of the same length as the `allow_types` list. `0` is a valid entry in the `num_overlap` list, as would be required if one wished to only specify overlap counts for a subset of the element types given in `allow_types`.

### Example 5c - Partial Overlap
In the same way that SVs can be made to completely overlap known element intervals in the various ways described above, they can also be made to *partially* overlap known element intervals. This can be done using SV config field `num_partial_overlap` which operated in the same way as `num_overlap`. An analogous version of the above example with partial overlap is given below:
```yaml
overlap_events:
    bed: ['/{path_to}/{candidate_overlap_events_1}.bed','/{path_to}/{candidate_overlap_events_2}.bed']
    allow_types: ['L1HS', 'L1PA3']
SVs:
    - type: "DEL"
      number: 10
      min_length: 500
      max_length: 1000
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
SVs:
    - type: "DIVERGENCE"
      number: 3
      divergence_prob: 0.2
      min_length:
        - 500
      max_length:
        - 1000
```
