# Output Description

## Output VCF file
The simulated variants are described in a VCF in [VCF 4.2 format](https://samtools.github.io/hts-specs/VCFv4.2.pdf).  
Predefined variants represented by a source region and a potential insertion locus are described 
by a single record.
More complex predefined variants: \
dDUP\_INV, INS\_iDEL dDUP\_iDEL, INVdel, delINV, INVdup, dupINV, dupINVdel, rTRA, INV\_rTRA, dupINVdel, delINVdup,
delINVdel, dupINVdup \
and Custom defined variants will be represented by a combination of INV, DEL and IDENTITY
operations.
Each operation can be in-place, meaning that the source region is modified or, pasted, in which case an insertion target has
to be provided to indicate where the potentially modified source region sequence has to be inserted.

VCF records relating to the same SV are connected by having the same value in the PARENT\_SVID INFO field, 
the PARENT\_SVTYPE INFO field gives the type of the original SV.
The original SV grammar is described in the GRAMMAR INFO field while the letter a record is affecting
appears in the SOURCE\_LETTER.

The details of the output representation for each built-in SV type are illustrated below.

## SVs represented as one VCF record

### INS: '' -> A

POS gives the reference base before which the novel sequence is inserted.
END is set equal to POS.   The inserted sequence is given in the INFO field INSSEQ.
The length of the inserted sequence is given as SVLEN.
The operation is marked paste in the IN\_PLACE field to indicate the affected letter does not 
appear in the left-hand side of the grammar.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	2	sv0	N	<INS>	100	PASS	END=2;SVTYPE=INS;VSET=0;IN_PLACE=paste;INSSEQ=GTTCGATGTACTCCCAGGCCGCTCATTCTCTCCAGTACTCAAAAGCAAGCTTTGC;SVLEN=55;INSORD=0;GRAMMAR=>A;SOURCE_LETTER=A 	GT	1|1
```

### DEL: A -> ''

POS gives the first deleted base.  END gives the last deleted base.  SVLEN gives the number of bases deleted.
The operation is marked in\_place in the IN\_PLACE field to indicate the source region was 
directly modified.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	6	sv0	N	<DEL>	100	PASS	END=60;SVTYPE=DEL;VSET=0;IN_PLACE=in_place;SVLEN=55;GRAMMAR=A>;SOURCE_LETTER=A 	GT	1|1
```

### INV: A -> a

POS gives the first base of the inverted region.  END gives the last base of the inverted
region.  SVLEN gives the length of the inverted region.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	6	sv0	N	<INV>	100	PASS	END=60;SVTYPE=INV;VSET=0;IN_PLACE=in_place;SVLEN=55;GRAMMAR=A>a;SOURCE_LETTER=A 	GT	1|1
```

### DUP: A -> AA

POS gives the first base of the duplicated region.  END gives the last base of the duplicated
region.  SVLEN gives the length of the duplicated region.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	6	sv0	N	<INV>	100	PASS	END=60;SVTYPE=INV;VSET=0;IN_PLACE=in_place;SVLEN=55;GRAMMAR=A>AA;SOURCE_LETTER=A 	GT	1|1
```
## SVs represented as multiple VCF records

### rTRA: A_B -> B_A
A reciprocal translocation is represented by a combination of four records.
A is pasted at B's starting base, this is represented by the first record of type IDENTITY, the insertion target
is given by the TARGET and TARGET\_CHROM fields. A is also deleted from its location in the reference,
this operation is represented by the first record of type DEL.
The two other records represent the symmetrical operations affecting B.
All the records share the same GRAMMAR, PARENT_SVID and PARENT_SVTYPE.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1	sv0_0	N	<IDENTITY>	100	PASS	END=21;SVTYPE=IDENTITY;VSET=0;IN_PLACE=paste;SVLEN=20;GRAMMAR=A_B>B_A;TARGET_CHROM=chr1;TARGET=45;PARENT_SVID=sv0;PARENT_SVTYPE=rTRA;SOURCE_LETTER=A 	GT	1|1
chr1	1	sv0_1	N	<DEL>	100	PASS	END=21;SVTYPE=DEL;VSET=0;IN_PLACE=in_place;SVLEN=20;GRAMMAR=A_B>B_A;PARENT_SVID=sv0;PARENT_SVTYPE=rTRA;SOURCE_LETTER=A 	GT	1|1
chr1	45	sv0_2	N	<IDENTITY>	100	PASS	END=60;SVTYPE=IDENTITY;VSET=0;IN_PLACE=paste;SVLEN=15;GRAMMAR=A_B>B_A;TARGET_CHROM=chr1;TARGET=1;PARENT_SVID=sv0;PARENT_SVTYPE=rTRA;SOURCE_LETTER=B 	GT	1|1
chr1	45	sv0_0	N	<DEL>	100	PASS	END=60;SVTYPE=DEL;VSET=0;IN_PLACE=in_place;SVLEN=15;GRAMMAR=A_B>B_A;PARENT_SVID=sv0;PARENT_SVTYPE=rTRA;SOURCE_LETTER=B 	GT	1|1
```

## Output PAF file
The correct whole-genome alignment of the simulated sample haplotypes to the reference can optionally be
written out in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).
To enable this, add the setting `output_paf: True`.   The alignment will
be written to files named `sim.hapA.paf` and `sim.hapB.paf`, representing the alignment
to the reference of `sim.hapA.fa` and `sim.hapB.fa`, respectively.











