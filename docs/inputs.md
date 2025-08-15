# Input Description
insilicoSV takes in a yaml configuration file as input, which specifies the path to the input reference genome and information regarding the SVs that will be simulated. 
Following the simulation, it outputs two haplotype files, and a VCF file describing the placement of the simulated SVs, a stats file, and a .fasta file containing any novel insertion sequences that were included in the simulation. 

### Parameter File
The configuration yaml file specifies simulation-wide settings (such as the reference genome) and a list of
variant sets to simulate.  For each variant set, it specifies the variant type, the range of sizes for each variant part,
and the number of variants to simulate. All configurations for variant sets should be put under the top-level "variant_sets" key. 
For example, the following configuration file defines three variant sets:

```yaml
# YAML config file
reference: {path}/{to}/ref.fa
variant_sets:
    - type: "INS"
      number: 10
      length_ranges: [[5, 10]]
    - type: "INVdel"
      number: 2
      length_ranges: 
        - [5, 10]
        - [5, 10]]
    - import: {path}/{to}/variants.vcf 
```

Variant sets can specify randomly simulated variants, or variants imported from a vcf.

The following parameters can be given for each randomly simulated variant set.  Each parameter is required unless otherwise specified.
1. *type*: str - the variant type or its corresponding grammar.  insilicoSV supports a predefined list of SV types and allows users to enter a custom transformation. Either a [valid grammar](sv_grammar.md) or one of the 26 predefined variant types named in [this figure](sv_grammar.md) should be entered.
2. *number*: int - describes how many of the specified variant type to simulate for this variant set.

For variants other than those related to tandem repeats (`trEXP` and `trCON`), the following
parameters must be given:
3. *length_ranges*: list of breakend distance ranges  - provides the minimum and maximum length for each variant part.  
SNPs, indels and simple SVs have a single part, while complex SVs
may have multiple parts; e.g. INVdel variants have two.  The order of part lengths in the list must correspond to the order of appearance
of part names used when defining the variant type; see [Example 2](use_cases.md#example-2---custom-svs) for an illustration. 
The values used for the ranges can be integers indicating a length in number of base pairs, or an expression relative to 
other parts of the SV (e.g., for an `rTRA`, `A_B -> B_A`, we can have `[[500, 1000], [0.5*A, 1.5*A]]` to indicate that interval 
B must be of length comprised between half and 1.5 times the length of A).
For predefined types with a dispersion, the length range of the dispersion will be in last position.
For Custom types, the length ranges are in order of appearance of the letters and the dispersions.
4. *overlap_mode [optional]*: str - enforce the SV to overlap a region defined in the files provided in `overlap_regions`. Must be `partial`, `contained`, `containing` or `exact` (see [example config](use_cases.md#example-5---placing-svs-into-specific-regions-of-interest-rois)).
5. *overlap_region_type [optional]*: list of str - only if an overlap mode is specified. Characterizes the regions to overlap, the name of the region has to contain one of the strings of the list. 

6. *interchromosomal: False [for SVs containing dispersions]*: Enable interchromosomal SVs. If True, each dispersion in the SV will be 
between two different chromosomes. All dispersions must be unbounded i.e. the dispersion range must be [null, null].
7. *interchromosomal_period: 0*: A value of `interchromosomal_period: 0` means each dispersion jumps to a different, randomly selected chromosome, while a value greater than 0 creates a cycle where the SV returns to the same set of chromosomes every `interchromosomal_period` dispersions.
8. *n_copies: [] [for SVs containing '+' grammar notation]*: specifies the number of copies for each sequence affected by a '+' in order of appearance in the grammar.
Each element of the list can be a positive number or a range of positive numbers. If a range is provided, a random number of copies included in the range will be used.
The default number of copies for a DUP is [1] and does not need to be specified.

For tandem repeat variants, the following parameters is needed:
8. *repeat_count_change_range*: the range from which to sample the number of repeats added (trEXP) or removed (trEXP).

For trEXP and trCON variants, a BED file of existing repeats must be specified in the 
*overlap_regions* global setting, and *overlap_region_type* for the existing repeat regions must
be specified in the variant set. The BED file columns must contain in order the chromosome, the start position, the end position, the region type, and the motif of the repeat.

SNPs and INDELs can be allowed to be overlapped by other SVs with the parameter:
9. *enable_overlap_sv [optional]*: bool - whether the SVs of this variant set can be overlapped by other SVs (default: False). Setting it to True
for SVs other than SNPs and INDELs will raise an error.

For specifying arm gain/loss or aneuploidy, the following parameters are available:
10. *arm_gain_loss=False [optional]*: bool - set to `True` for the SV to duplicate or delete an entire chromosome arm.
11. *arm_percent=[100, 100] [optional]*: [int, int] - a range [min_percent, max_percent] from which to determine the percentage of the chromosome arm to duplicate or delete. The operation starts from the extremity of the arm. 
12. *aneuploidy=False [optional]*: bool - set to True for the SV to duplicate or delete a whole chromosome copy.
13. *aneuploidy_chrom=None [optional]*: List[str] - A list of chromosome names (e.g., ['chr21', 'chrX']) on which aneuploidy is permitted to occur. If None, aneuploidy can occur on any chromosome.

The following parameters can be set on the top level of the config file and provide higher-order controls over SV placement:
1. *reference*: str - path to input reference used as template for simulation.
2. *blacklist_regions*: list[str] - list of paths to BED or VCF files containing intervals to be ignored during SV placement (see [example config](use_cases#example-4---marking-banned-intervals-of-the-genome)).
3. *overlap_regions*: list[str] - list of paths to BED files containing genome elements to be used for overlapping SV placement (see [example config](use_cases#example-5---placing-svs-at-known-repetitive-element-intervals)).
4. *max_tries: 100 [optional]*: int - number of tries to find a valid position to simulate each SV.
5. *max_random_breakend_tries: 100 [optional]*: int - number of tries to find a breakend by taking a random position in the genome before checking available regions.
6. *homozygous_only: False [optional]*: bool - if set to True, make all simulated variants homozygous
7. *heterozygous_only: False [optional]*: bool - if set to True, make all simulated variants heterozygous
8. *filter_small_chr [optional]*: int - filter out chromosomes of length less than the given integer (if no value is provided then no filtering will occur).
9. *th_proportion_N: 0.05 [optional]*: The proportion of N and n base pairs an SV is allowed to cover.
10. *enable_hap_overlap=False [optional]*: Enable heterozygous SVs to overlap across homologous chromosomes.
11. *arms [optional]*: str - A path to a file specifying the centromere positions for the chromosomes. Refer to the  [use_cases](use_cases.md#example-7---chromosome-gainloss)
for a detailed description of the required file format.

Examples of the full set of simulation options available through various config inputs can be found in the [use cases](use_cases) page.

### Importing SVs from VCF files
insilicoSV offers the possibility to edit the reference with known SVs.
An example of the syntax is shown in the config file above and [this example](use_cases.md#example-3---editing-reference-with-known-svs)

insilicoSV currently supports importation of any VCF in insilicoSV [VCF format](outputs.md#output-vcf-file) or 
VCF files in [Gnomad](https://gnomad.broadinstitute.org/) format containing simple SVs only, meaning, all the ones
not reported in the CPX_TYPE field and not reported as BND.
Other types would have to be adapted to insilicoSV format [VCF format](outputs.md#output-vcf-file).

When writing a VCF in insilico format for SV importation, the required fields are the CHROM, POS and ID fields as well as
the INFO field with an END, SVTYPE, SVLEN, and, if applicable, TARGET, TARGET_CHROM, SVID. 

### BED files format
The BED files must contain the following columns in order: Chromosome, Start Position,
End Position, Region Name.
BED files for tandem repeat regions must additionally contain a fifth column with the Motif of the region.
Additional columns can appear after those required columns.
