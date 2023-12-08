# Input and Output Description
insilicoSV takes in a yaml configuration file as input, which specifies the path to the input reference genome and information regarding the SVs that will be simulated. Following the simulation, it outputs two haplotype files, BED and VCF files describing the placement of the simulated SVs, a stats file, and a .fasta file containing any novel insertion sequences that were included in the simulation. 

### Parameter File
The configuration yaml file specifies the path to the input reference genome as well as the SV type, range of lengths, and the number of SVs of each type to simulate. All configurations for SVs should be put under the "variant_sets" key (note: "SV" here includes SNPs and INDELs because while they are not technically types of SVs, they are treated in the same way by the simulator). For each SV, the following parameters are available:
1. *type*: str - insilicoSV supports a predefined list of SVs and allows users to enter a custom transformation. Either "Custom" or one of the 19 predefined SV types given in [SV grammar](https://github.com/PopicLab/insilicoSV-dev/blob/develop/docs/sv_grammar.md) should be entered.
2. *number*: int - describes how many of the specified SV type to simulate
3. *min_length*: int - If an SV only has a single reference interval then a single integer should be given for the min_length and max_length. If an SV has multiple reference intervals then multiple must be provided to min_length and max_length to specify the size ranges for each component of the SV. Each min_length / max_length entry is allocated to the SV reference intervals in lexicographical order (see [Example 2](https://github.com/PopicLab/insilicoSV-dev/blob/develop/docs/example_use_cases.md#example-2---custom-svs) for an illustration).
4. *max_length*: int - must be the same type as min_length, note that max_length >= min_length >= 0 for all elements in each
5. *source=None [optional]* - Source sequence for a custom SV, see below to find instructions on how to create a transformation
6. *target=None [optional]* - Target sequence for a custom SV, see below to find instructions on how to create a transformation

The following parameters can be set under the "sim_settings" key to change default configurations for the simulator:
1. *reference*: str - path to input reference used as template for simulation
2. *max_tries=100 [optional]*: int - number of tries to find a valid position to simulate each SV
3. *fail_if_placement_issues=False [optional]*: bool - if set to True, insilicoSV will raise an Exception when a single SV fails to be placed and simulated
4. *filter_small_chr [optional]*: int - filter out chromosomes of length less than the given integer (if no value is provided then no filtering will occur).
5. *prioritize_top=False [optional]*: bool - if set to True, variants will be added to the reference in the order in which they are listed in the config. This prioritizes the variants listed first because if the reference becomes overcrowded and there is no more remaining space to place non-overlapping SVs, then the remaining SVs will not be included in the output. If set to False, the full set of SVs will be shuffled and added to the reference in random order.
6. *homozygous_only=False [optional]*: bool - if set to True, make all simulated variants homozygous

The following parameters can be set on the top level of the config file and provide higher-order controls over SV placement:
1. *avoid_intervals*: str - path to VCF containing intervals to be ignored during SV placement (see [example config](https://github.com/PopicLab/insilicoSV-dev/blob/develop/docs/example_use_cases.md#example-4---marking-banned-intervals-of-the-genome))
2. *overlap_events*: str - path to BED file containing genome elements to be used for overlapping SV placement (see [example config](https://github.com/PopicLab/insilicoSV-dev/blob/develop/docs/example_use_cases.md#example-5---placing-events-at-known-repetitive-element-intervals))

A basic example configuration file is given below, and examples of the full set of simulation options available through various config inputs can be found in the [example use cases](https://github.com/PopicLab/insilicoSV-dev/blob/develop/docs/example_use_cases.md) page.
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
SVs:
    - type: "DEL"
      number: 3
      min_length: [1000]
      max_length: [10000]
```

### Output BED File
Each line/entry of the BED file describes a single SV component, which we describe as an *event*, meant to indicate the fundamental mutation operators that compose in different ways to constitute SVs of different type. Each BED line will have the following parameters:
1. *source_chr*: The source chromosome of the event
2. *source_start*: Start position on the source_chr
3. *source_end*: End position on the source_chr
4. *target_chr*: The target chromosome of the event [*note: currently only intra-chromosomal SVs are supported*]
5. *target_start*: Start position on the target chr
6. *target_end*: End position on the target chr
7. *event_type*: Describes the transformation made by the event, either an INS, DEL, INV, TRA, DUP, INVDUP, or INVTRA.
8. *event_size*: Size of the reference fragment impacted by the event
9. *zygosity*: {0/1, 1/0} = heterozygous, 1/1 = homozygous. insilicoSV gives each SV a 50% chance of being heterozygous or homozygous, and if the SV is heterozygous it is given a 50% chance of being placed on haplotype A or B. 
10. *parent_sv*: Describes the parent SV the event is a component of (for instance, the INV record from a dupINVdup will have "dupINVdup" as its parent SV). If a custom SV was provided, the name becomes "source>target"
11. *nth_sv*: int, index to count each SV. All events of a SV belong in the same index.
12. *order*: int, for insertion-like operations such as TRA, INS, or DUP, the "order" index describes in which order events that target the same position were compiled. Events with INV and DEL operations have an order of 0.

Example BED file output are given below:
```
chr1    109334938   109335727   chr1    109338930   109338930   DUP 789 1/1 dDUP    1   1
chr7    130007589   130007849   chr7    130007589   130007849   DEL 260 1/1 DEL 2254    0
chr12	85434890    85435083	chr12	85434890    85435083	INV 193	0/1 INV	3751	0
```

### Output VCF file
The output BED file is accompanied by a VCF describing the same set of SVs but in the [VCF 4.2 format](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Whereas the BED file entries each describe individual sub-SV events, the VCF entries each describe an entire SV. One augmentation made to the output VCF format is the use of an info field called `TARGET` to describe the target locus of  a dispersion-based SV. For instance, the VCF record for a dispersed duplication is given with start and end values that describe the SV's source interval (i.e., the interval that is duplicated), and the `TARGET` field records the position at which the copy is inserted into the reference.

Example VCF outputs (of the same SVs given in the BED file output) are given below:
```
chr1    109334939   dDUP    N   <dDUP>  100 PASS    END=109335727;SVTYPE=dDUP;SVLEN=789;TARGET=109338930    GT  1/1
chr7    130007589   130007849   chr7    130007589   130007849   DEL 260 1/1 DEL 2254    0
chr12   85434891    INV N   <INV>   100 PASS    END=85435083;SVTYPE=INV;SVLEN=193   GT  0/1
```
