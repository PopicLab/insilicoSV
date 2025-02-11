# SV Grammar

insilicoSV uses a grammatical notation to represent how SVs of a given type transform the reference genome.
This notation can be used to simulate custom SV types involving arbitrary transformations.
It is also used to give precise definitions for built-in SVs, as illustrated below:

![Graphical illustration of insilicoSV grammar](gallery/main.png)

An SV type is represented by a single grammar rule, where the left-hand side (LHS) represents one or more
reference intervals, while the right-hand side (RHS) represents their rearrangement and/or transformation.  In
the LHS, capital letters represent reference intervals; letters adjacent in the LHS represent adjacent intervals
in the reference.  Each letter on the LHS must appear there exactly once.
In the RHS, a capital letter represents the sequence of the corresponding reference
interval, while its lowercase version represents the reverse complement of that sequence.
Letters present in the LHS but not in the
RHS represent deletions, while letters present in the RHS but not in the LHS represent novel insertions.
Thus, A $\rightarrow$ a represents a simple inversion, while AB $\rightarrow$ b represents an inversion
flanked by a deletion.

Each RHS letter present in the LHS represents either an in-place transform of a reference interval,
or an insertion of the interval's sequence at a new location. Thus, AB $\rightarrow$ BA
represents a swap of two adjacent intervals, effectuated by translocating the first interval immediately
downstream of the second.   

Duplications of a reference interval are represented by including
the corresponding letter more than once in the RHS.  For example, A $\rightarrow$ AA represents
a simple duplication, while AB $\rightarrow$ Aba represents a duplication-flanked inversion.

Adding a + after a letter on the right-hand side represents a variable number
of copies of the corresponding reference interval: for example, A $\rightarrow$ AA+ represents adding one or more
copies of A.

Dispersions can be represented using underscores.  Thus, a downstream dispersed duplication can be written
as A_ $\rightarrow$ A_A, and a reciprocal translocation of two intervals as A_B $\rightarrow$ B_A.

An illustration of using the grammar to specify custom SVs is given in
[this example](use_cases#example-2---custom-svs).