# SV Grammar

insilicoSV uses a grammatical notation to represent how SVs of a given type transform the reference genome.
This notation can be used to simulate custom SV types involving arbitrary transformations.
It is also used to give precise definitions for built-in SVs, as illustrated below:

![Graphical illustration of insilicoSV grammar](sample_imgs/fig1a.pdf)

An SV type is represented by a single grammar rule, where the left-hand side (LHS) represents one or more
reference regions, while the right-hand side (RHS) represents their rearrangement and/or transformation.  In
the LHS, capital letters represent reference regions; letters adjacent in the LHS represent regions adjacent
in the reference.  Each letter on the LHS must appear there exactly once.
In the RHS, a capital letter represents the sequence of the corresponding reference
region, while its lowercase version represents the reverse complement of that sequence.
Letters present in the LHS but not in the
RHS represent deletions, while letters present in the RHS but not in the LHS represent novel insertions.
Thus, A $\rightarrow$ a represents a simple inversion, while AB $\rightarrow$ b represents an inversion
flanked by a deletion.

Each RHS letter present in the LHS represents either an in-place transform of a reference region,
or an insertion of the region's sequence at a new location.  Those RHS letters representing
an insertion at a new location are marked by appending a caret (^).  Thus, AB $\rightarrow BA^
represents a swap of two adjacent regions, effectuated by translocating the first region immediately
downstream of the second.   The caret markings serve to disambiguate among alternate operational
representations of an SV.  For example, AB $\rightarrow BA^ denotes the same transformation of the
reference, but represented as an upstream rather than a downstream translocation.
For RHS letters denoting novel sequence, the caret mark is always implied.

Duplications and inverted duplications of a reference region are represented by including
the corresponding letter more than once in the RHS.  For example, A $\rightarrow$ AA^ represents
a simple duplication, while AB $\rightarrow$ Aba^ represents a duplication-flanked inversion.
Note that at most one instance of a given letter on the RHS may lack a caret mark.

Adding a + after a letter on the right-hand side represents a variable number
of copies of the corresponding reference region: for example, A -> AA+ represents adding one or more
copies of A.  For RHS letters with a +, the caret mark is always implied.

Dispersions can be represented using underscores.  Thus, a downstream dispersed duplication can be written
as A_ -> A_A^, and a reciprocal translocation of two regions as A_B -> B^_A^.

An illustration of using the grammar to specify custom SVs is given in
[this example](example_use_cases.md#example-2---custom-svs).
