# Input Description
```insilicoSV``` takes in a YAML configuration file as input, which specifies the path to the input reference genome and information regarding what SVs to simulate. 

### Parameter File
The configuration YAML file specifies global parameters (e.g., the reference genome) and the definition of one or multiple variant categories to simulate. The definition of each variant category/set includes the variant type, a list of breakend distance ranges (specifying the minimum and maximum distance allowed between SV breakends), variant placement contraints,
and the number of such variants to simulate. All variant definitions should be specified under the top-level "variant_sets" parameter. Variant set definitions can specify both randomly simulated variants or variants imported from a VCF. For example, the following configuration file defines three variant categories:

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

The following parameters can be specified for each randomly-simulated variant set.  Each parameter is required unless otherwise specified.
1. *type*: str - the variant type or its corresponding grammar.  insilicoSV supports a predefined list of SV types and allows users to enter a custom transformation. Either a [valid grammar](sv_grammar.md) or one of the 26 predefined variant types named in [this figure](sv_grammar.md) should be entered.
2. *number*: int - describes how many of the specified variant type to simulate for this variant set.

For variants other than tandem repeats (`trEXP` and `trCON`), the following parameters must be given:

3. *length_ranges*: list of breakend distance ranges  - provides the minimum and maximum reference interval length allowed between breakends.  
SNPs, INDELs and simple SVs consist of a single reference interval, while complex SVs may involve multiple intervals; e.g. INVdel variants have two. The order of the interval lengths in the list must correspond to the order of appearance of their corresponding letter in the variant type definition; see [Example 2](use_cases.md#example-2---custom-svs) for an illustration. The values used to define the range can be integers indicating the number of base pairs, or an expression relative to other SV intervals (e.g., for an `rTRA`, `A_B -> B_A`, we can have `[[500, 1000], [0.5*A, 1.5*A]]` to indicate that interval B must be of length comprised between half and 1.5 times the length of A).
For predefined types with a dispersion, the length range of the dispersion will be in the last position.
For Custom types, the length ranges are in order of appearance of the letters and the dispersions.

4. *overlap_mode [optional]*: str - enforce the SV to overlap a region defined in the files provided in `overlap_regions`. Must be `partial`, `contained`, `containing`, `exact`, `terminal` or `whole-chromosome` (see [example config](use_cases.md#example-5---placing-svs-into-specific-regions-of-interest-rois)).

5. *overlap_region_type [optional]*: list of str - only if an overlap mode is specified. Characterizes the regions to overlap, the name of the region has to contain one of the strings of the list. 

6. *interchromosomal: False [for SVs containing dispersions]*: Enable interchromosomal SVs. If True, each dispersion in the SV will be 
between two different chromosomes. All dispersions must be unbounded i.e. the dispersion range must be [null, null].
7. *interchromosomal_period: 0*: A value of `interchromosomal_period: 0` means each dispersion jumps to a different, randomly selected chromosome, while a value greater than 0 creates a cycle where the SV returns to the same set of chromosomes every `interchromosomal_period` dispersions.
8. *n_copies: [] [for SVs containing '+' grammar notation]*: specifies the number of copies for each sequence affected by a '+' in order of appearance in the grammar.
Each element of the list can be a positive number or a range of positive numbers. If a range is provided, a random number of copies included in the range will be used.
The default number of copies for a DUP is [1] and does not need to be specified.

For tandem repeat variants, the following parameters is needed:

9. *repeat_count_change_range*: the range from which to sample the number of repeats added (trEXP) or removed (trEXP).

For trEXP and trCON variants, a BED file of existing repeats must be specified in the 
*overlap_regions* global setting, and *overlap_region_type* for the existing repeat regions must
be specified in the variant set. The BED file columns must contain in order the chromosome, the start position, the end position, the region type, and the motif of the repeat.

SNPs and INDELs can be placed within SVs using the parameter:

10. *allow_sv_overlap [optional]*: bool - set to `True` to allow this variant set to be overlapped by SVs (default: False). Setting it to True
for variants other than SNPs and INDELs will raise an error.

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
10. *allow_hap_overlap=False [optional]*: Allow heterozygous SVs to overlap across homologous chromosomes.

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
