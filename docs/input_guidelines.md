# Input Guidelines
insilicoSV takes in two input files: the reference genome and a yaml configuration file. Following the simulation, it outputs two haplotype files, a BED file, and a simple stats file using the prefix given. 

### Reference Genome
This file should be in FASTA format.

### Parameter File
The configuration yaml file specifies the range of the lengths of the SVs along with the number to simulate. All configurations for structural variants should be put under the "SVs" key. For each SV, the following parameters are available (with some additional parameters available for certain specific use cases described in [the next section](example_use_cases.md)):
1. *type*: str, insilicoSV supports a predefined list of SVs and allows users to enter a custom transformation. Either "Custom" or one of the 16 predefined SV types named in the below table should be entered.
2. *number*: int, describes how many of the specified SV type to simulate
3. *min_length*: int or list, if an integer is provided, insilcoSV assumes that each event's length within a SV must fall between the min_length and max_length. Entering a list offers customization by specifying a different range for each event. If provided a list, insilicoSV assumes lengths are entered to correspond with the symbols in lexicographical order. The lengths for non-dispersion events (represented by alphabetical characters) will therefore be considered before that of dispersions (represented by an underscore). See [usage example #2](example_use_cases.md#example-2---custom-svs). 
4. *max_length*: int or list, must be the same type as min_length, note that max_length >= min_length >= 0 for all elements in each
5. *source=None [optional]*: Source sequence for a custom SV, see below to find instructions on how to create a transformation
6. *target=None [optional]*: Target sequence for a custom SV, see below to find instructions on how to create a transformation

The following optional parameters can be set under the "sim_settings" key to change default configurations for the simulator:
1. *max_tries=100 [optional]*: number of tries to find a valid position to simulate each SV
2. *fail_if_placement_issues=False [optional]*: if set to True, insilicoSV will raise an Exception when a single SV fails to be placed and simulated
3. *generate_log_file=False [optional]*: if set to True, insilicoSV will generate a log file for diagnostic purposes and debugging
4. *filter_small_chr={int} [optional]*: filter out chromosomes of length less than the given integer
5. *prioritize_top=False [optional]*: if set to True, simulate the SVs listed first as the later ones may fail to be placed. Otherwise, give each SV an equal chance to be simulated
6. *homozygous_only=False [optional]*: if set to True, make all simulated variants homozygous


### Output BED File

Each line/entry will have the following parameters:
1. *source_chr*: The source chromosome of the event
2. *source_start*: Start position on the source_chr [INCLUDE at pos], zero-based indexing
3. *source_end*: End position on the source_chr [EXCLUDE at pos], one-based indexing
4. *target_chr*: The target chromosome of the event
5. *target_start*: Start position on the target chr [INCLUDE at pos], zero-based indexing
6. *target_end*: End position on the target chr [EXCLUDE at pos], one-based indexing
7. *event_type*: Describes the transformation made by the event, either an INS, DEL, INV, TRA, DUP, INVDUP, or INVTRA. Dispersed duplications--those that do not occur immediately after the original--have an attached "d" at the front.
8. *event_size*: Size of the reference fragment impacted by the event
9. *zygosity*: {0/1, 1/0} = heterozygous, 1/1 = homozygous. insilicoSV gives each SV a 50% chance of being heterozygous or homozygous, and if the SV is heterozygous it is given a 50% chance of being placed on haplotype A or B. 
10. *parent_sv*: Describes the parent SV the event is a component of, for instance "dupINVdup." If a custom SV was provided, the name becomes "source>target"
11. *nth_sv*: int, index to count up each SV (note: not the events). All events of a SV belong in the same index.
12. *order*: int, for insertion-like operations such as TRA, INS, or DUP, the "order" index describes in which order events that target the same position were compiled. Events with INV and DEL operations have an order of 0.

Example BED file output are given below:
```
chr1    109334938   109335727   chr1    109338930   109338930   DUP 789 1/1 dDUP    1   1
chr7    130007589   130007849   chr7    130007589   130007849   DEL 260 1/1 DEL 2254    0
chr12	85434890    85435083	chr12	85434890    85435083	INV 193	0/1 INV	3751	0
```

### Output VCF file
The output BED file is accompanied by a VCF describing the same set of SVs but in the (VCF 4.2 specification)[https://samtools.github.io/hts-specs/VCFv4.2.pdf]. One augmentation made to the output VCF format is the use of an info field called `TARGET` to describe the target locus of  a dispersion-based event. For instance, the VCF record for a dispersed duplication is given with start and end values that describe the SV's source interval (i.e., the interval that is duplicated), and the `TARGET` field records the position at which the copy is inserted into the reference.

Example VCF outputs (of the same events given in the BED file output) are given below:
```
chr1    109334939   dDUP    N   <dDUP>  100 PASS    END=109335727;SVTYPE=dDUP;SVLEN=789;TARGET=109338930    GT  1/1
chr7    130007589   130007849   chr7    130007589   130007849   DEL 260 1/1 DEL 2254    0
chr12   85434891    INV N   <INV>   100 PASS    END=85435083;SVTYPE=INV;SVLEN=193   GT  0/1
```

**Note**: In the case of an overcrowded simulation (i.e., so many variants specified for too little total reference sequence length) not all SVs inputted may be simulated. Because the simulator randomly selects positions for SVs, some SVs may fail to be placed due to oversaturation. 
