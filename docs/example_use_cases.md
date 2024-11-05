# Example Use Cases

insilicoSV provides various simulation features that can be used together or separately to generate
synthetic genomes with varying levels of control over SV placement. Examples of the different use cases with
matching config file structure are provided below. For all of these cases, commandline usage is:

```
insilicosv -c <config.yaml>
```

### Example 1 - Predefined SV types

To incorporate SVs from the built-in library of types, a configuration file of the following form can be
provided with parameters provided for the count and min/max length for each set of SVs to be included in the
output genome.

```yaml
# YAML config file
sim_settings:
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
    - type: "dupINVdel"  # "ABC" -> "Aba^"
      number: 1
      length_ranges:
        - [5, 10]  # min/max length for first DEL
        - [10, 15]  # min/max length for INV
        - [5, 10]  # min/max length for second DEL
```

### *Example 1a* - Example config with entire insilicoSV vocabulary
This [summary config](summary_config.md) contains an example specifying the event size ranges for each predefined variant type,
as well as examples of the various placement features that can be used to bias variants towards or away from regions of interest.

### Example 1b - Example SNP/INDEL specification
SNPs and INDELs are supported by insilicoSV as well. Because SNPs are only a single base in length, they only need to be specified with `number` as in the example below:
```yaml
variant_sets:
  - type: "SNP"  # "A" -> "A*" (for A of length 1)
    number: 10
```
INDELs are not given a unique type label but can be simulated by setting a sufficiently small min and max size for an SV of type DEL or INS.

### Example 1c - Unbounded dispersions
When specifying length ranges for dispersion intervals, `None` can be provided as an upper bound and this will result in
the dispersion target locus being placed at any position on the same chromosome as the source event. For example:`
```yaml
variant_sets:
  - type: "dDUP"  # "A_" -> "A_A^" or "_A" -> "A^_A"
    number: 2
    length_ranges:
      - [500, 1000]  # min/max length for source
      - [500, None]  # unbounded dispersion: target locus may be placed arbitrarily far
```
Unbounded dispersions can be used with predefined or custom SVs.

### Example 2 - Custom SVs
Custom SVs can be specified by manually describing the desired variant with the grammatical notation described in [SV grammar](sv_grammar.md). An example input config is given below:
```yaml
# YAML config file
sim_settings:
    reference: "{path}/{to}/ref.fa"
variant_sets:
    - type: "Custom"
      source: "AB_C_D"
      target: "bb^_A^Ec_EDC^"
      number: 1
      length_ranges: 
        - [5, 10]   # A
        - [6, 10]   # B
        - [7, 10]   # C
        - [8, 10]   # D
        - [10, 15]  # E
        - [10, 15]  # first _
        - [15, 20]  # second _
```

In the "target", caret markings indicate "source" symbol instances which represent insertions of (possibly inverted)
sequence of the source symbol at a new location, rather than an in-place transformation of the reference region represented
by the source symbol.  At most one instance of each source symbol in the target may be left without a caret.  The order of
caret-less symbols in the target must match their order in the source.

The source and target must have the same numbers of dispersions, with the i'th dispersion in the source
corresponding to the i'th dispersion in the target.  As one exception, all dispersions in the source
can be omitted, in which case  symbols in the target are used to infer the target section
corrsponding to the source.  For example, with `source: "ABC"`, `target: "A\_B^\_C^"`
causes source to be interpreted as `source: "ABC__"`, while `target: "A^\_B\_C^"` causes
source to be interpreted as `source: "_ABC_"`.

### Example 3 - Editing reference with known SVs
To edit an input reference file with a known set of variants the user can provide them in a VCF file. 
The predefined variant types supported for this use case are DEL, DUP, INV, INS, dDUP, INV\_dDUP, TRA\_UNBALANCED, INV\_DUP3,
and SNP. For insertions, VCF records may be specified with the insertion sequence given in an INFO field called `INSSEQ`
(provided a matching header line is included as well). All VCF records are expected to include an info field `SVTYPE`
to record variant type. The config file to perform this reference edit is:
```yaml
# YAML config file
sim_settings:
    reference: "{path}/{to}/ref.fa"
variant_sets:
    - import: "{path_to_vcf}"
```

If a population-level VCF is being provided in which records include an allele frequency info field (see
[VCF spec. 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf)), each variant's genotype will be sampled
based on the allele frequency; variants for which the sampling yields the reference allele for both
haplotypes will be omitted from import.

Variant sets imported from vcf may be combined in the same config file with variant sets specifying
random variants to simulate.

### Example 4 - Marking genome blacklist intervals
When initializing a new simulation the user can include a blacklist of genome intervals (i.e., intervals that specified
variants will avoid) via a VCF or BED file given (or multiple provided in a comma-separated list) in the `blacklist_regions`
entry of the config file:
```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
blacklist_regions: "{path}/{to}/{banned_intervals}.{bed, vcf}"
variant_sets:
    - type: "DEL"  # "A" -> ""
      number: 3
      length_ranges: [[1000, 10000]]
      blacklist_region_type: "all"
    ...
```
Each variant set specified in the config can include a `blacklist_region_type` value which will control which intervals
recorded in `blacklist_regions` will be avoided by those variants. In the above example, `"all"` is specified, which will
result in none of the three deletions from being placed in any of the regions in `blacklist_regions`. If `blacklist_regions`
is given input in BED file format, records can include a fourth column recording `type`, which can then be used to filter
the blacklist intervals considered for a given set of SVs.

If blacklist entries are provided in VCF format, only the following parts of the record are used for the blacklist: CHROM, POS and
the INFO fields END, SVLEN, TARGET and TARGET_CHROM.

### Example 4b - Specifying a minimum inter-variant distance
Variant placement can also be constrained by enforcing that there be a minimum inter-variant distance between any two
breakends belonging to different variants. This value will be by default 1 (i.e., enforcing that any two variants
are separated by at least 1 bp) but a desired simulation-specific distance can be provided in the config input `min_intersv_dist`:
```yaml
sim_settings:
  reference: "{path}/{to}/ref.fa"
  min_intersv_dist: 100  # <- each DEL will be separated by at least 100bp
variant_sets:
  - type: "DEL"
    number: 3
    length_ranges: [(1000, 10000)]
```


### Example 5 - Placing SVs at known regions of interest

To augment a randomized simulation of SVs onto an input reference, the user can include in the
simulation config file the path to a BED file (or multiple) containing known element intervals (e.g.,
known repetitive elements taken from RepeatMasker).  If the `overlap_mode` field is given a value
(`"exact"`, `"partial"` or `"contained"`), then each variant of that set will be placed in exact or
partial overlap with a randomly selected interval from `overlap_regions`.  An example config with
these inputs is:

```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
overlap_regions: ["/{path_to}/{candidate_overlap_events}.bed"]
variant_sets:
    - type: "DEL"  # "A" -> ""
      number: 5
      length_ranges: [[null, null]]
      overlap_mode: "exact"  # <- or "partial" or "contained"
```

In `"exact"` overlap mode, the SV must exactly match the overlap region on which it is placed.
The length of the SV is then determined by the overlap region, and must be left unspecified
in `"length_ranges"`.  In `"partial"` overlap mode, the SV must overlap the overlap region
partially (but not completely).  In `"contained"` overlap mode, the SV must be fully contained
within the overlap region.

Multiple BED files can be given as input and their records will be combined and drawn from during SV placement (in this
case the user should provide a list of paths). Records from the `overlap_regions` BED file(s) will be shuffled on input,
and the file is required to have the first four columns of standard BED records (chrom, chromStart, chromEnd, name).  Files
specifying known Tandem Repeat regions for expansion/contraction need to have a fifth column specifying the repeat unit
length of each repeat region.

Regions relevant to placing variants from a given variant set can be filtered down by region type
and length.  `overlap\_region\_type` specifies a list of prefixes; regions whose name (4th bed
column) starts with one of the prefixes will be used.  `overlap\_region\_length\_range` can be
specified to give a range (a 2-element [min, max] list) of regions to include.  Either min or max can
be `null` to leave that side of the range open.  (Note that specifying overlap region length range is
entirely separate from specifying length ranges for SV components.)

The output VCF file will label which SVs were placed at specified intervals with the additional INFO field
`OVLP={evt. name}', as in this example record:
```
chr21   18870078    DEL N   DEL 100 PASS    END=18876908;SVTYPE=DEL;SVLEN=6831;OVLP=L1HS  GT  0/1
```

### Example 5a - Placing specific SV components at regions of interest

The portion of the SV which participates in overlap with an overlap region is termed an _anchor_.
For SVs without dispersions, the anchor defaults to the full SV.  To constrain the placement of SVs with
dispersions, or to specify an anchor other than the full SV, the anchor can be specified as part
of the SV's grammar definition.  For example:

```yaml
sim_settings:
    reference: "{path}/{to}/ref.fa"
overlap_regions: ["/{path_to}/{candidate_overlap_events_1}.bed"]
variant_sets:
    - type: "delINVdel"  # "ABC" -> "b"
      source: "A(B)C"
      number: 5
      length_ranges:
        - [500, 1000]
        - [null, null]
        - [500, 1000]
      overlap_mode: "exact"
      overlap_region_type: ["L1HS"]
```

Parentheses indicate the part(s) of the SV that are constrained to overlap with an overlap region according
to the overlap mode.  The anchor can be placed on either the source or the target, and can wrap
any contiguous sub-sequence of source or target elements (including an empty one),
subject to the following limitations: (1) the anchor cannot contain dispersions;
(2) for `"exact"` overlap mode, the anchor must contain exactly one SV component, whose length range must
be left unspecified.  An empty anchor is only compatible with the overlap mode `"contained"`.

### Example 6 - Divergent intervals
In addition to the various SV types included in the predefined library, insilicoSV also allows for the simulation of divergent
intervals in which a random proportion of the bases in a given interval will be corrupted with some probability `p`. Divergent
intervals can be included in the same way as any of the SVs provided, accessible through the built-in type name `'DIVERGENCE'` (or through custom SVs using the asterisk divergence operator as described in [SV grammar](sv_grammar.md)). The
probability parameter by which each base will be randomly changed can be optionally provided by the user as shown in the example
below:
```yaml
sim_settings:
    reference: {path}/{to}/ref.fa
variant_sets:
    - type: "DIVERGENCE"  # "A" -> "A*"
      number: 3
      divergence_prob: 0.2
      length_ranges: [[500, 1000]]
```
