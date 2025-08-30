# Use Cases

insilicoSV provides various simulation features that can be used together or separately to generate
synthetic genomes with varying levels of control over SV simulation and placement. Examples of the different 
use cases with matching config files are provided below. 

### Example 1 - Predefined SV types

To incorporate SVs from the built-in library of types, a configuration file of the following form can be
provided with parameters provided for the count and length ranges for each set of SVs to be included in the
output genome. Note that the length of the dispersion must always be the last entry in length_ranges.

```yaml
# YAML config file
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
    - type: "rTRA"  # "A_B" -> "B_A"
      number: 1
      length_ranges:
        - [5, 10]  # A
        - [10, 15]  # NOTE: B
        - [5, 10]  # NOTE: _
```

### *Example 1a* - Example config with entire insilicoSV vocabulary
This [summary config](summary_config.md) contains an example specifying the event size ranges for each predefined variant type,
as well as examples of the various placement features that can be used to bias variants towards or away from regions of interest.

### Example 1b - Example SNP/INDEL specification
SNPs and INDELs are supported by insilicoSV as well. Because SNPs are only a single base in length, they only need to be specified with `number` as in the example below:
```yaml
variant_sets:
  - type: "SNP" 
    number: 10
```
An `INDEL` can be defined using the pre-defined INDEL type, or by specifying DEL or INS types with a length range below 50.
If you choose the pre-defined `INDEL` type, each variant within the set will be randomly assigned as either a DEL or an INS. 
You can also omit the `length_ranges` parameter for the INDEL type, in which case it will automatically default to a range of `[[1, 50]]`.

### Example 1c - Unbounded dispersions
When specifying length ranges for dispersion intervals, `null` can be provided as an upper bound and this will result in
the dispersion target locus being placed at any position on the same chromosome as the source event. For example:`
```yaml
variant_sets:
  - type: "dDUP"  # "A_" -> "A_A" or "_A" -> "A_A"
    number: 2
    length_ranges:
      - [500, 1000]  # min/max length for source
      - [500, null]  # unbounded dispersion: target locus may be placed arbitrarily far
```
Unbounded dispersions can be used with predefined or custom SVs.

### Example 2 - Custom SVs
Custom SVs can be specified by manually describing the desired variant with the grammatical notation described 
in [SV grammar](sv_grammar.md). Note that for custom SVs the length_ranges of the letters AND the dispersions must be provided
in their order of appearance from left to right.
Note that the same number of dispersions must appear in the lhs and the rhs.
An example input config is given below:
```yaml
# YAML config file
reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "AB_C_D->bb_AEc_EDC"
      number: 1
      length_ranges: 
        - [5, 10]   # A
        - [6, 10]   # B
        - [10, 15]  # first _
        - [7, 10]   # C
        - [15, 20]  # second _
        - [8, 10]   # D
        - [10, 15]  # E
```

The source and target must have the same numbers of dispersions, with the i'th dispersion in the source
corresponding to the i'th dispersion in the target.  

N.B. the grammar is not symmetrical, for instance, in the example below:
```yaml
variant_sets:
    - type: "A_B->B_A"
      number: 1
      length_ranges: 
        - [5, 10]   # A
        - [10, 15]  # _
        - [15, 20]   # B
```
The interval B, of size in the range [15, 20], will be after the interval A in the reference, of size in the range [5, 10].
Thus, this custom event is not symmetrical.
However, we defined the predefined SV types as symmetrical. Thus, for instance, rTRA corresponds to
"A_B" -> "B_A" or "B_A" -> "A_B" the order being selected randomly. The previous example defined as a rTRA would 
allow the smaller interval A to be after the larger interval B.

### Example 3 - Editing reference with known SVs
To edit an input reference file with a known set of variants the user can provide them in a VCF file. 
The predefined variant types supported for this use case are DEL, DUP, INV, INS, dDUP, INV\_dDUP, rTRA, 
DUP\_INV, and SNP. For insertions, VCF records may be specified with the insertion sequence given in an INFO 
field called `INSSEQ`
(provided a matching header line is included as well). All VCF records are expected to include an info field `SVTYPE`
to record variant type. The config file to perform this reference edit is:
```yaml
# YAML config file
reference: "{path}/{to}/ref.fa"
variant_sets:
    - import: "{path_to_vcf}"
```

The homozygosity or heterozygosity of SVs is determined using the GT field. A genotype of 1|1 or 1/1 will be considered
as homozygous while 0/1, 0|1, 1/0, and 1|0 will be considered as heterozygous and the modified haplotype will be 
randomly selected.
If a sample with a 'GT' field is not specified, the genotype will be chosen at random.

Variant sets imported from VCFs may be combined in the same config file with variant sets specifying
random variants to simulate.

### Example 4 - Marking genome blacklist intervals
When initializing a new simulation the user can include a blacklist of genome intervals 
(i.e., intervals that specified variants will avoid) via a VCF or BED file given (or multiple provided 
in a comma-separated list) in the `blacklist_regions`
entry of the config file:
```yaml
reference: "{path}/{to}/ref.fa"
blacklist_regions: "{path}/{to}/{banned_intervals}.{bed, vcf}"
variant_sets:
    - type: "DEL"  # "A" -> ""
      number: 3
      length_ranges: [[1000, 10000]]
      blacklist_region_type: "all"
    - type: "INS"  # "" -> "A"
      number: 3
      length_ranges: [[100, 1000]]
```
Each variant set specified in the config can include a `blacklist_region_type` value which will control which intervals
recorded in `blacklist_regions` will be avoided by those variants. In the above example, `"all"` is specified for the DELs, 
which will result in none of the three deletions from being placed in any of the regions in `blacklist_regions`. 
If `blacklist_regions` is given in BED file format, records can include a fourth column recording `type`, 
which can then be used to filter the blacklist intervals considered for a given set of SVs.
If `blacklist_regions` is given in VCF file format, records can include the info field `REGION_TYPE` which can then be used to filter the blacklist intervals considered for a given set of SVs.
If `REGION_TYPE` is not provided, the name of the region will be `DEFAULT`.
In this example, the three insertions will be placed randomly regardless of the blacklist regions as the blacklist_region_type
is not specified.

If blacklist entries are provided in VCF format, only the following parts of the record are used for the blacklist: 
CHROM, POS and the END INFO field.

### Example 4b - Specifying a minimum inter-variant distance
Variant placement can also be constrained by enforcing that there be a minimum inter-variant distance between any two
breakends belonging to different variants. This value will be by default 1 (i.e., enforcing that any two variants
are separated by at least 1 bp) but a desired simulation-specific distance can be provided in the config input `min_intersv_dist`:
```yaml
reference: "{path}/{to}/ref.fa"
min_intersv_dist: 100  # <- any two SVs will be separated by at least 100bp
variant_sets:
  - type: "DEL"
    number: 3
    length_ranges: [[1000, 10000]]
  - type: "INS"
    number: 3
    length_ranges: [[100, 1000]]
```


### Example 5a - Placing SVs into specific regions of interest (ROIs)

To constrain the placement of SVs to specific regions, the path to a single or multiple BED files containing these intervals (e.g.,
known repetitive elements taken from RepeatMasker) can be provided.  Setting the `overlap_mode` field in a specific 
variant set to either `"exact"`, `"partial"`, `"containing"` or `"contained"`, will enforce that each variant of that set 
will be placed (with the corresponding overlap mode) into a randomly selected interval from `overlap_regions`.  
An example config with these inputs is:

```yaml
reference: "{path}/{to}/ref.fa"
overlap_regions: ["/{path_to}/{candidate_overlap_events}.bed"]
variant_sets:
    - type: "DEL"  # "(A)" -> ""
      number: 5
      length_ranges: [[null, null]]
      overlap_mode: "exact"  # <- or "partial" or "contained"
```

In `"exact"` overlap mode, the SV must exactly match the ROI on which it is placed.
The length of the SV is then determined by the ROI, and must be left unspecified
in `"length_ranges"`.  In `"partial"` overlap mode, the SV must overlap the ROI
partially (but not completely). In `"containing"` overlap mode, the SV must strictly contain
the ROI. In `"contained"` overlap mode, the SV must be strictly contained
within the ROI.

Multiple BED files (provided as a list of paths) can be given as input and their records will be combined and drawn 
from during SV placement. Each file is required to have the first four columns of standard BED records 
(chrom, chromStart, chromEnd, name).  Files specifying known Tandem Repeat regions for expansion/contraction 
need to have a fifth column specifying the motif of each repeat region.

ROIs relevant to placing variants from a given variant set can be filtered down by region type
and length.  `overlap\_region\_type` specifies a list of identifiers; regions whose name (4th bed
column) containing one of the identifiers will be used.  `overlap\_region\_length\_range` can be
specified to give overlap constraints (a 2-element [min, max] list).  Either min or max can
be `null` to leave that side of the range open.  (Note that specifying an overlap length range is
entirely separate from specifying length ranges for SV components.)

The output VCF file will label which SVs were placed at specified intervals with the additional INFO field
`OVLP={region name}', as in this example record:
```
chr21   18870078    DEL N   DEL 100 PASS    END=18876908;SVTYPE=DEL;SVLEN=6831;OVLP=L1HS;VSET=0;IN_PLACE=in_place;GRAMMAR=A>AA;SOURCE_LETTER=A  GT  0/1
```

For SVs without dispersions, the anchor defaults to the whole SV.  To constrain the placement of SVs with
dispersions, or to specify an anchor other than the full SV, the anchor can be specified as part
of the SV's source grammar definition.  For example:

```yaml
reference: "{path}/{to}/ref.fa"
overlap_regions: ["/{path_to}/{candidate_overlap_events_1}.bed"]
variant_sets:
    - type: "A(BC) -> b"  # delINVdel
      number: 5
      length_ranges:
        - [500, 1000]
        - [null, null]
        - [null, null]
      overlap_mode: "exact"
      overlap_region_type: ["L1HS"]
```

Parentheses indicate which part(s) of the SV are constrained to overlap with an ROI according
to the overlap mode.  The anchor must be placed on the source, and can wrap
any contiguous sub-sequence of source elements (including an empty one).
For `"exact"` overlap mode, the length ranges of the constrained SV's part(s) must
be left unspecified.  An empty anchor is only compatible with the overlap mode `"contained"` and
can be used to constrain an insertion target for instance.

If an anchor constrains a dispersion, the dispersion will be intrachromosomal even if the SV is set as interchromosomal 
(In this case, other unconstrained dispersions of the SV will then be interchromosomal). 


### Example 5b - Inversions within segmental duplications (SDs)
INVs that span the regions between SDs can be simulated using the `"exact"` overlap placement mode given a BED file containing these regions:
```yaml
reference: "{path}/{to}/ref.fa"
overlap_regions: ["/{path_to}/{regions_between_SDs}.bed"]
variant_sets:
    - type: "INV"
      number: 5
      length_ranges: [[null, null]]
      overlap_mode: "exact"
```

To obtain the required BED file, `bedtools complement` can be used to generate regions that are not covered by SDs (note: depending on the use case, the terminal regions should be removed as a post-processing step):
```
bedtools complement -i your_sd_regions.bed -g your_reference/ref.fa > regions_between_SDs.bed
```

### Example 6 - Interchromosomal Dispersions and Interchromosomal Periods
You can use the `interchromosomal` and `interchromosomal_period` flags to control how an SV is dispersed across chromosomes. 
By default, `interchromosomal` is False, which means the SV is intrachromosomal (it stays on the same chromosome).

#### Understanding the `interchromosomal_period` Flag
The `interchromosomal_period` flag gives you fine-grained control over which chromosomes are affected by an interchromosomal SV.

- `interchromosomal_period=0`: This is the default.  Each dispersion of the SV will switch chromosome. 
This ensures that every part of the SV is placed on a different chromosome than the last.

- `interchromosomal_period > 0`: This creates a cycle of dispersions among a specific number of different chromosomes.
The first `interchromosomal_period` dispersions will each be placed on a new, unique chromosome.
Any subsequent dispersions will then cycle back through these same chromosomes in the order they were first visited creating a cycle.

- List of Integers (e.g., [2, 3]): You can provide a list of two integers to define a range. 
For each SV in the variant set, a random `interchromosomal_period` value will be chosen from this range.

##### Key Considerations for Interchromosomal Dispersions
- `Unbounded Lengths`: When defining interchromosomal dispersions, their lengths must be unbounded (specified as [null, null]).  
- Simplified Syntax: If you provide a non-null `interchromosomal_period` value in your configuration, the SV will automatically be treated as `interchromosomal`, 
even if the `interchromosomal` flag is not explicitly set to True.
- If `interchromosomal_period=0`, a chromosome can be repeated, however a dispersion is ensured to connect two different chromosomes.

#### Example scenario: Cycling between Chromosomes
The provided YAML code defines a variant set that places cycles of templated insertions across chromosomes.
```yaml
reference: "{path}/{to}/ref.fa"
overlap_regions: ["/{path_to}/{candidate_overlap_events_1}.bed"]
variant_sets:
    - type: "A_B_C -> A_BCAB_C" 
      number: 5
      length_ranges:
        - [500, 1000]
        - [null, null]
        - [500, 1000]
        - [null, null]
        - [500, 1000]
      interchromsomal: True
      interchromosomal_period: [2, 3]
```
In this example, for each of the 5 structural variants defined, the `interchromosomal_period` will be randomly set to either 2 or 3. 
This means some SVs will cycle between two chromosomes, while others will cycle between three, creating a diverse set of interchromosomal insertions.

### Example 7 - Chromosome Gain/Loss
- Unbounded Lengths: When defining interchromosomal dispersions, their lengths must be unbounded (specified as [null, null]).  
- Multiple Dispersions: If a custom SV with multiple dispersions is flagged as interchromosomal, each dispersion will involve a change to a different chromosome.
- Example Scenario: In the second example provided above, where the type is "A__ -> A_AB_A", a possible placement for this interchromosomal SV could be:
  - The source of A is on chr1.
  - The AB segment is placed on chr3.
  - The second copy of A is placed on chr2.
  
  Because all dispersions are interchromosomal, the last copy of A cannot be placed back on chr3. 
  However, it could be placed in a different region of chr1.

### Example 7 - SNP and INDEL placement within SVs
SNPs and INDELs can be allowed to overlap SVs as follows:

```yaml
reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "DUP"  
      number: 5
      length_ranges:
        - [1000, 100000]
    - type: "INDEL"
      number: 10
      length_ranges:
         - [30, 50]
      allow_sv_overlap: True
    - type: "SNP"
      number: 20
      allow_sv_overlap: True
```
Here the INDEL and SNP definitions have the `allow_sv_overlap` parameter set to `True`, which allows them to be randomly placed within the DUP intervals.
Note: SNPs and INDELs that are allowed to overlap SVs are always considered as occurring first in the simulation process. 
As such, they might be modified or even deleted by SVs that are placed later. Regardless of whether they are ultimately observable in the final genome, all simulated variants are included in the final VCF output.

### Example 8 - Chromosome Gain/Loss
This section details parameters for simulating chromosome arm gain/loss or whole chromosome aneuploidy.

#### Arm Gain/Loss (`arm_gain_loss: True`)
To enable the duplication or deletion of entire chromosome arms, set the `arm_gain_loss` parameter to `True`.

##### Centromere File (arms)
- **Purpose:** A BED file containing the centromere start and end positions for each chromosome of your reference genome. This file is required when arm_gain_loss is True.
- **File Format:** The first four columns of the BED file must be:
  - Chromosome name 
  - Chromosome length 
  - Beginning of the centromere (start coordinate)
  - End of the centromere (end coordinate)
- **Scope:** The file does not need to include all reference chromosomes; only those for which arm gain/loss is intended.

**arm_percent parameter:** specifies a range (e.g., `[60, 80]`) 
to determine the percentage of the chromosome arm to be duplicated or deleted, starting from the arm's extremity.

#### Aneuploidy (aneuploidy: True)
**aneuploid_chrom parameter:** You may optionally provide a list of specific chromosome names (e.g., `['chr1', 'chr22']`) 
on which aneuploidy is permitted. If this parameter is not provided, aneuploidy can occur on any chromosome in the reference.

#### General Considerations
- **Compatible SV Types:** Only DEL (Deletion) and DUP (Duplication) are compatible with arm_gain_loss and aneuploidy flags.
- **length_ranges:** When simulating aneuploidy, length_ranges should either not be provided or be set to `[[null, null]]`.
- **Multiple Aneuploid DUPs:** Multiple aneuploid duplications can target the same chromosomes to increase chromosome copy numbers arbitrarily. In such cases, overlap constraints will be disregarded.
- **Heterozygous Nature:** Chromosome gain/loss and aneuploidy events are defined as heterozygous, affecting only one of the existing chromosome copies.
- **New Chromosome Copies (DUP Aneuploidy):** For a DUP with `aneuploidy: True`, new chromosome copies will be created and named chrom_copy_num (where num is the copy number).
  - If n_copies is not provided or set to 1, a case of trisomy (one additional chromosome copy) will be simulated.
  - As shown in the example, `n_copies: 3` allows for the creation of three additional chromosome copies.

```yaml
reference: "{path}/{to}/ref.fa"
arms: ["/{path_to}/{arm_regions}.bed"] # Required for arm_gain_loss, but not for aneuploidy
variant_sets:
    - type: "DUP" 
      number: 3
      aneuploidy: True
      aneuploid_chrom: ['chr1', 'chr22']
      n_copies: 3 # Simulates 3 additional copies for selected chromosomes
    - type: "DEL" 
      number: 5
      arm_gain_loss: True
      arm_percent: [60, 80] # Deletes 60-80% of a chromosome arm from its extremity
```
