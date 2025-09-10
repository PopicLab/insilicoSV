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
If the pre-defined `INDEL` type is chosen, each variant within the set will be randomly assigned as either a DEL or an INS. 
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

### Example 4a - Marking genome blacklist intervals
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
You can constrain where an SV is placed by providing a list of genomic intervals in one or more **BED files** under the 
`overlap_regions` parameter.  You can then specify an `overlap_mode` for a given variant set to control how each SV interacts with these regions.

### Overlap modes and constraints
There are several ways to define the relationship between an SV and a region of interest (ROI):
* **`"exact"`**: The constrained breakends of the SV have to exactly match the ROI boundaries. The SV's length is therefore 
determined by the ROI, so `[null, null]` must be used for the corresponding range in `length_ranges`.
* **`"partial"`**: The constrained breakends of the SV must overlap with one of the boundaries of a selected ROI.
* **`"containing"`**: The constrained breakends of the SV must completely contain a selected ROI.
* **`"contained"`**: The constrained breakends of the SV must be completely contained within a selected ROI.
* **`"terminal"`**: The constrained breakends of the SV is placed at the extremity of a chromosome arm.
* **`"whole-chromosome"`**: The SV spans an entire chromosome. This mode is only compatible with **Deletions (DEL)** and **Duplications (DUP)**. 
For Duplications, setting `n_copies` to `1` (or not specifying it) creates a single additional chromosome copy (trisomy).

### Defining and using anchors
An **anchor** defined by parentheses `()` in the SV's source grammar can be used to specify which part of a complex SV should overlap with an ROI.

For example, in `A(BC) -> b`, the `BC` portion is the anchor. It is the only part of the SV that will be constrained to 
overlap with a given ROI. The anchor can be any contiguous sub-sequence of the source. For SVs without dispersions, 
the entire SV defaults to being the anchor.

An **empty anchor** `()` can be used to constrain the target location of an insertion and is only compatible with the `"contained"` overlap mode.


### General considerations

* **Filtering Regions**: To filter the ROIs from BED files based on their name (4th column) using `overlap_region_type` 
or by length using `overlap_region_length_range`.
* **N-Regions**: By default, insilicoSV avoids placing SVs in regions with a high proportion of 'N' base pairs (over 10%). 
This can affect placements in telomeres (`overlap_mode: terminal`). You can adjust this threshold using the global parameter `th_proportion_N`.
* **Intrachromosomal Constraint**: If a dispersion is part of an anchor, it will be forced to be intrachromosomal, even if the overall SV is marked as interchromosomal.
The potential remaining dispersions will be interchromosomal.
* **BED Format**: Multiple BED files (provided as a list of paths) can be given as input and their records will be combined and drawn 
from during SV placement. Each file is required to have the first four columns of standard BED records 
(chrom, chromStart, chromEnd, name).  Files specifying known Tandem Repeat regions for expansion/contraction 
need to have a fifth column specifying the motif of each repeat region.
* **Output VCF**: The final VCF file will include an `OVLP` field in the `INFO` column for each SV placed within a specified region, indicating the name of the region.

#### Example YAML configuration

```yaml
reference: "{path}/{to}/ref.fa"
overlap_regions: ["/{path_to}/{candidate_overlap_events_1}.bed"]
variant_sets:
    - type: "AB(C) -> b"  # delINVdel
      number: 5
      length_ranges:
        - [500, 1000]
        - [500, 1000] 
        - [null, null] # Inside the anchor of an `exact` overlap
      overlap_mode: "exact"
      overlap_region_type: ["L1HS"] # Allowed regions for overlap
```
In this example, 5 SVs of type `delINVdel` are created. Their `C` region is constrained to overlap with a region from 
the `candidate_overlap_events_1.bed` file, with the overlap mode set to `"exact"`.

#### Example output VCF file
```
chr21   18870078    DEL N   DEL 100 PASS    END=18876908;SVTYPE=DEL;SVLEN=6831;OVLP=L1HS;VSET=0;IN_PLACE=in_place;GRAMMAR=A>;SOURCE_LETTER=A  GT  0/1
```
#### Simulating chromosome arm-level SVs
Overlapping constraints can be used to simulate SVs affecting an entire chromosome arm. 
This is particularly useful for modeling large-scale genetic events like whole-arm deletions or duplications.
For this purpose, the precise start and end coordinates for each chromosome arm have to be provided in a BED file. 
Then, applying the overlap_mode `exact` to a variant set will enforce the resulting SVs to completely
span a chromosome arm.

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
bedtools complement -i sd_regions.bed -g reference/ref.fa > regions_between_SDs.bed
```

### Example 6 - Interchromosomal dispersions and interchromosomal periods
The `interchromosomal` and `interchromosomal_period` parameters control how an SV is dispersed across chromosomes. 
By default, `interchromosomal` is False, which means the SV is intrachromosomal (it stays on the same chromosome).

#### Understanding the `interchromosomal_period` flag
The `interchromosomal_period` flag gives you fine-grained control over which chromosomes are affected by an interchromosomal SV.

- `interchromosomal_period=0`: This is the default if `interchromosomal` is set to True.  Each dispersion of the SV will switch chromosome. 
This ensures that every part of the SV is placed on a different chromosome than the last.

- `interchromosomal_period > 0`: This creates a cycle of dispersions among a specific number of different chromosomes.
The first `interchromosomal_period` dispersions will each be placed on a new, unique chromosome.
Any subsequent dispersions will then cycle back through these same chromosomes in the order they were first visited creating a cycle.

- List of Integers (e.g., [2, 3]): A list of two integers can be provided to define a range. 
For each SV in the variant set, a random `interchromosomal_period` value will be chosen from this range.

##### Key considerations for interchromosomal dispersions
- Unbounded lengths: When defining interchromosomal dispersions, their lengths must be unbounded (specified as [null, null]).  
- Simplified syntax: If a non-null `interchromosomal_period` value is provided, the SV will automatically be treated as `interchromosomal`, 
even if the `interchromosomal` flag is not explicitly set to True.
- If `interchromosomal_period=0`, a chromosome can be repeated, however a dispersion is ensured to connect two different chromosomes.

#### Example scenario: cycling between chromosomes
The provided YAML code defines a variant set that places cycles of templated insertions across chromosomes.
```yaml
reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "A_B_C -> A_BCAB_C" 
      number: 5
      length_ranges:
        - [500, 1000]
        - [null, null]
        - [500, 1000]
        - [null, null]
        - [500, 1000]
      interchromosomal: True
      interchromosomal_period: 1
```
In this example, for each of the 5 structural variants defined, the `interchromosomal_period` will be set to 1. 
This means these SVs will cycle between two chromosomes.

For instance, a valid placement for one of these SVs would be:
- A is on chr2.
- B is on chr1.
- C is on chr2.

The VCF records for this case would be:
```
chr1   74348760    sv0_2   N   <COPY-PASTE>    100 PASS    END=74349707;OP_TYPE=COPY-PASTE;GRAMMAR=A_B_C->A_BCAB_C;VSET=0;TARGET_CHROM=chr1;TARGET=74349707;SVLEN=947;INSORD=2;SVID=sv0;SVTYPE=CUSTOM;SYMBOL=B    GT  0|1
chr2   49674622    sv0_1   N   <COPY-PASTE>    100 PASS    END=49675453;OP_TYPE=COPY-PASTE;GRAMMAR=A_B_C->A_BCAB_C;VSET=0;TARGET_CHROM=chr1;TARGET=74349707;SVLEN=831;INSORD=1;SVID=sv0;SVTYPE=CUSTOM;SYMBOL=A    GT  0|1
chr2   101202093   sv0_0   N   <COPY-PASTE>    100 PASS    END=101202831;OP_TYPE=COPY-PASTE;GRAMMAR=A_B_C->A_BCAB_C;VSET=0;TARGET_CHROM=chr1;TARGET=74349707;SVLEN=738;INSORD=0;SVID=sv0;SVTYPE=CUSTOM;SYMBOL=C   GT  0|1
```

Another valid example without period:
```yaml
reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "A_B_C -> A_BCAB_C" 
      number: 5
      length_ranges:
        - [500, 1000]
        - [null, null]
        - [500, 1000]
        - [null, null]
        - [500, 1000]
      interchromosomal: True
```
In this case, `interchromosomal` has been set to `True` without providing any `interchromosomal_period`.
Therefore, the period will default to 0 and each jump will swicth chromosomes without any cycle constraints.
Note: a chromosome might be visited several times.

A possible placement for one of these interchromosomal SVs could be:
  - A is on chr1.
  - B is on chr13.
  - C is on chr20.
  
Because all dispersions are interchromosomal, B cannot be placed on chr1 and C cannot be placed on chr13. 
However, C could be placed in a different region of chr1.

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

### Example 8 - Divergence
The divergence `*` symbol can be used to introduce point mutations in a duplicated sequence, causing it to differ from the original reference sequence. 
The rate of mutations can be configured using the `divergence_prob` parameter, which specifies the probability of each nucleotide in the sequence to undergo a point mutation.

```yaml
reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "A->AA*"  
      number: 5
      divergence_prob: 0.1
      length_ranges:
        - [1000, 1000]
```
The provided example defines a divergent tandem DUP of 1kbp length. Each of the 1kbp duplicated nucleotides will have a `10%` chance of being mutated.
