# SV Grammar

insilicoSV uses a grammatical notation to represent the various types of SVs that can be supported by the simulator. A graphical representation of the notation is given below:

[[sample_imgs/fig.png]]

The form of a given SV's transformation of the genome is represented with a grammatical notation that describes the mapping from input reference sequence to output synthetic donor sequence. With this notation, the reference sequence is given by a series of capital letters corresponding to the intervals affected by the SV, and the donor sequence is given by a transformation of those letters that reflects the ways in which the reference intervals are mutated by the SV. For example, a deletion-flanked inversion (delINV) is notated as AB $&#8594$ b, in which the left side of the expression indicates the two reference intervals involved in the delINV and the right side indicates the donor sequence that will appear in place of AB, that is the inverted interval b (the lowercase indicating that B will appear inverted in the donor). Although SNPs are not considered to be a type of structural variant, we include them here as another valid event type for simulation (see [use cases](example_use_cases.md#example-1b---example-snp-specification) for usage examples).

A custom SV consists of a user-generated transformation with a source and target sequence of symbols. Some examples of a source would be ABC and A_B_C, while some examples of the target would be a'AB or A_b'Bc'_C. 

insilicoSV maps a random fragment of the reference to each of the symbols in the source and recompiles the affected region with the target. For instance, AB $&#8594$ A would remove the fragment marked as symbol B. All symbols in the source sequence MUST be unique to create a one-to-one mapping between symbol and reference fragment.

| Name | Symbol | Description |
|------|--------|-------------|
| Generic Event | Any uppercase alphabetical letter | The most fundamental organizing tool that maps to a reference fragment |
| Inversion | Any lowercase alphabetical letter | Indicates an inversion. <br /> Ex. a transformation ABC $&#8594$ abc will invert A, B, and C and organize the new fragments as denoted in the target |
| Duplication | Original symbol followed by single quotation (') | An original symbol refers to the initial character used in the source sequence. *There can only be ONE original symbol for every unique character - all other copies, including those that are inverted, must have a duplication marking (').* <br /> Ex. A transformation ABC $&#8594$ ABA'c would duplicate A after B and invert the fragment C.|
| Dispersion | Underscore (_) | Indicates a gap between the symbols surrounding it. Note that events may be simulated within a dispersion but not within other events. |
| Insertions | Uppercase alphabetical letter | To add foreign, randomly-generated insertions, use a symbol not present in the source to the target sequence. <br /> Ex. A_B $&#8594$ A_BC inserts a randomly-generated sequence after the fragment indicated by symbol B|
