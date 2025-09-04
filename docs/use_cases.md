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


### Example 5 - Placing SVs into specific regions of interest (ROIs)
You can constrain where an SV is placed by providing a list of genomic intervals in one or more **BED files** under the 
`overlap_regions` parameter.  You can then specify an `overlap_mode` for a given variant set to control how each SV interacts with these regions.

### Overlap Modes and Constraints
There are several ways to define the relationship between an SV and a region of interest (ROI):

* **`"exact"`**: The constrained breakends of the SV have to exactly match the ROI boundaries. The SV's length is therefore 
determined by the ROI, so you must use `[null, null]` for the corresponding range in `length_ranges`.
* **`"partial"`**: The constrained breakends of the SV must overlap with one of the boundaries of a selected ROI.
* **`"containing"`**: The constrained breakends of the SV must completely contain a selected ROI.
* **`"contained"`**: The constrained breakends of the SV must be completely contained within a selected ROI.
* **`"terminal"`**: The constrained breakends of the SV is placed at the extremity of a chromosome arm.
* **`"whole-chromosome"`**: The SV spans an entire chromosome. This mode is only compatible with **Deletions (DEL)** and **Duplications (DUP)**. 
For Duplications, setting `n_copies` to `1` (or not specifying it) creates a single additional chromosome copy (trisomy).

### Defining and Using Anchors

To specify which part of a complex SV should overlap with an ROI, you can use an **anchor** defined by parentheses `()` in the SV's source grammar.

For example, in `A(BC) -> b`, the `BC` portion is the anchor. It is the only part of the SV that will be constrained to 
overlap with a given ROI. The anchor can be any contiguous sub-sequence of the source. For SVs without dispersions, 
the entire SV defaults to being the anchor.

An **empty anchor** `()` can be used to constrain the target location of an insertion and is only compatible with the `"contained"` overlap mode.


### General Considerations

* **Filtering Regions**: You can filter the ROIs from your BED files based on their name (4th column) using `overlap_region_type` 
or by length using `overlap_region_length_range`.
* **N-Regions**: By default, InsilicoSV avoids placing SVs in regions with a high proportion of 'N' base pairs (over 10%). 
This can affect placements in telomeres (`overlap_mode: terminal`). You can adjust this threshold using the global parameter `th_proportion_N`.
* **Intrachromosomal Constraint**: If a dispersion is part of an anchor, it will be forced to be intrachromosomal, even if the overall SV is marked as interchromosomal.
The potential remaining dispersions will be interchromosomal.
* **BED Format**: Multiple BED files (provided as a list of paths) can be given as input and their records will be combined and drawn 
from during SV placement. Each file is required to have the first four columns of standard BED records 
(chrom, chromStart, chromEnd, name).  Files specifying known Tandem Repeat regions for expansion/contraction 
need to have a fifth column specifying the motif of each repeat region.
* **Output VCF**: The final VCF file will include an `OVLP` field in the `INFO` column for each SV placed within a specified region, indicating the name of the region.

#### Example YAML Configuration

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

#### Example Output VCF File
```
chr21   18870078    DEL N   DEL 100 PASS    END=18876908;SVTYPE=DEL;SVLEN=6831;OVLP=L1HS;VSET=0;IN_PLACE=in_place;GRAMMAR=A>;SOURCE_LETTER=A  GT  0/1
```

#### Simulating Chromosome Arm-Level SVs
You can use overlapping constraints to simulate SVs affecting an entire chromosome arm. 
This is particularly useful for modeling large-scale genetic events like whole-arm deletions or duplications.

For this purpose, you need a BED file that defines the precise start and end coordinates for each chromosome arm. 
Then, you can apply the overlap_mode `exact` to the variant sets of your choice to enforce the resulting SVs to completely
span a chromosome arm.

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

### Example 7 - Placing overlapping SVs
Any SNP or INDEL can be allowed to be overlapped by other, non-overlapping SVs.

```yaml
reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "nrTRA"  
      number: 5
      length_ranges:
        - [500, 1000]
        - [500, 1000]
    - type: "DEL"
      number: 10
      length_ranges:
         - [30, 50]
      enable_overlap_sv: True
    - type: "SNP"
      number: 20
      enable_overlap_sv: True
```
The DEL and SNP sets have the `enable_overlap_sv: True` setting, which allows them to be overlapped by the nrTRA variants.
DEL and SNP variants cannot overlap each other.

SVs that are allowed to overlap are always considered as occurring first in the simulation process. 
This means they might be modified or even deleted by other SVs that are placed later. 
Regardless of whether they are ultimately observable in the final genome, all overlapping SVs will be included in the final VCF output.

### Example 8 - Time Point Mutations for Cancer Genome Modeling
For more complex SV interactions, such as those found in cancer lineages, insilicoSV provides a workflow in the notebook 
`insilicosv_cancer_genome.ipynb`. This notebook is designed to model cancer evolution, generate genomes for different 
clonal cell populations, and simulate reads at varying tumor purity levels.

The notebook requires a configuration file that specifies the reference genome, read coverage, and paths to
`insilicoSV` configuration files. These configuration files describe the SVs at different time points, the clones 
corresponding to each time point, and the purity of each clone. Example configuration files are available in the 
`configs/clone_configs/` folder. 

For more detailed information, you should refer to the notebook itself.
