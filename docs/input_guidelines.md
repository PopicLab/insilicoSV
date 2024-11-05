# Input and Output Description
insilicoSV takes in a yaml configuration file as input, which specifies the path to the input reference genome and information regarding the SVs that will be simulated. Following the simulation, it outputs two haplotype files, BED and VCF files describing the placement of the simulated SVs, a stats file, and a .fasta file containing any novel insertion sequences that were included in the simulation. 

### Parameter File
The configuration yaml file specifies simulation-wide settings (such as the reference genome) and a list of
variant sets to simulate.  For each variant set, it specifies the variant type, the range of sizes for each variant part, and the number of variants to simulate. All configurations for variant sets should be put under the top-level "variant_sets" key. For example, the following configuration file defines two variant sets:

```yaml
# YAML config file
sim_settings:
    reference: {path}/{to}/ref.fa
variant_sets:
    - type: "INS"
      number: 10
      length_ranges: [[5, 10]]
    - type: "INVdel"
      number: 2
      length_ranges: [[5, 10], [5, 10]]
```

Variant sets can specify randomly simulated variants, or variants imported from a vcf.

The following parameters can be given for each randomly simulated variant set.  Each parameter is required unless otherwise specified.
1. *type*: str - the variant type.  insilicoSV supports a predefined list of SV types and allows users to enter a custom transformation. Either "Custom" or one of the 19 predefined variant types named in [this figure](sv_grammar.md) should be entered.
2. *number*: int - describes how many of the specified variant type to simulate for this variant set.

For variants other than those related to tandem repeats (`trINS`, `trEXP` and `trCON`), the following
parameters must be given:
3. *length_ranges*: list of 2-element lists - provides the minimum and maximum length for each variant part.  SNPs, indels and simple SVs have a single part, while complex SVs
may have multiple parts; e.g. INVdel variants have two.  The order of part lengths in the list must correspond to the alphabetical order
of part names used when defining the variant type; see [Example 2](example_use_cases.md#example-2---custom-svs) for an illustration. List with one 2-element range should be provided for variants with a single source interval.
5. *source=None [for variant sets of type Custom]*: Source sequence for a custom SV, see below to find instructions on how to create a transformation
6. *target=None [for variant sets of type Custom]*: Target sequence for a custom SV, see below to find instructions on how to create a transformation

For tandem repeat variants, the following parameters are needed:
7. *repeat_count_change_range*: the range from which to sample the number of repeats added or removed

For trINS variants, the following parameters are needed:
8. *repeat_sequence_choices*: list of strings from which the single repeat sequence will be chosen

For trEXP and trCON variants, a .bed file of existing repeats must be specified in the 
*overlap_regions* global setting, and *overlap_region_type* for the existing repeat regions must
be specified in the variant set.  The repeat unit length of each existing repeat needs to be
given in the 5th column of the .bed file.


The following parameters can be set under the "sim_settings" key to change default configurations for the simulator:
1. *reference*: str - path to input reference used as template for simulation
2. *max_tries=100 [optional]*: int - number of tries to find a valid position to simulate each SV
6. *homozygous_only=False [optional]*: bool - if set to True, make all simulated variants homozygous

The following parameters can be set on the top level of the config file and provide higher-order controls over SV placement:
1. *blacklist_regions*: list[str] - list of paths to BED or VCF files containing intervals to be ignored during SV placement (see [example config](example_use_cases.md#example-4---marking-banned-intervals-of-the-genome))
2. *overlap_regions*: list[str] - list of paths to BED files containing genome elements to be used for overlapping SV placement (see [example config](example_use_cases.md#example-5---placing-svs-at-known-repetitive-element-intervals))

Examples of the full set of simulation options available through various config inputs can be found in the [example use cases](example_use_cases.md) page.

### Output VCF file
The simulated variants are described in a VCF in [VCF 4.2 format](https://samtools.github.io/hts-specs/VCFv4.2.pdf).  Each variant is represented as a combination of basic SVs of one of these built-in types:
'DEL', 'INS', 'INV', 'DUP', 'INV\_DUP', 'INV\_DUP2', 'dDUP', 'INV\_dDUP', 'TRA\_NONRECIPROCAL', 'INV\_TRA'.
vcf records relating to the same SV are connected by having the same value in the PARENT\_SVID INFO field,
and the PARENT\_SVTYPE INFO field gives the type of the original SV.

### Output PAF file
The correct whole-genome alignment of the simulated sample haplotypes to the reference can optionally be
written out in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).
To enable this, add the setting `output_paf: True` under `sim_settings`.   The alignment will
be written to files named `sim.hapA.paf` and `sim.hapB.paf`, representing the alignment
to the reference of `sim.hapA.fa` and `sim.hapB.fa`, respectively.

