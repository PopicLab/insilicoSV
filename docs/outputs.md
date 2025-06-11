# Output Description

## Output VCF file
The simulated variants are described in a VCF in [VCF 4.2 format](https://samtools.github.io/hts-specs/VCFv4.2.pdf).  
Predefined variants represented by a source region and a potential insertion locus are described 
by a single record.
More complex predefined variants: \
dDUP\_INV, INS\_iDEL dDUP\_iDEL, INVdel, delINV, INVdup, dupINV, dupINVdel, rTRA, INV\_rTRA, dupINVdel, delINVdup,
delINVdel, dupINVdup \
and Custom defined variants will be represented by a combination of CUT, COPY-PASTE, CUT-PASTE, COPYinv_PASTE, CUTinv-PASTE, and INV 
operations.
Each operation can be cut-pasted, meaning that the source region is deleted or, copy-pasted, in which case the source region is conserved.
In both cases, an insertion target is provided to indicate where the potentially modified source region sequence has to be inserted.
An INV can modify a copy-paste or cut-paste operation to indicate that the pasted sequence is inverted or, INV can be an operation
on its own.
A CUT operation refers to a deletion of a sequence without it being pasted somewhere else.
The operation a VCF record corresponds to is indicated in the OP_TYPE INFO field.

VCF records relating to the same SV are connected by having the same value in the SVID INFO field, 
the SVTYPE INFO field gives the type of the SV.
The SV grammar is described in the GRAMMAR INFO field while the letter a record is affecting
appears in the SYMBOL INFO field.
In addition, for SVs represented by multiple records, the OP_TYPE field characterizes the type of 
transformation performed on the source.

The details of the output representation for each built-in SV type are illustrated below.

## SVs represented as one VCF record

### INS: '' -> A

POS gives the reference base before which the novel sequence is inserted.
END is set equal to POS.   The inserted sequence is given in the INFO field INSSEQ.
The length of the inserted sequence is given as SVLEN.
The OP_TYPE is NA indicating it is a simple SV and does not need to be decomposed in different operations.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	2	sv0	N	<INS>	100	PASS	END=2;OP_TYPE=NA;GRAMMAR=->A;VSET=0;INSSEQ=GTTCGATGTACTCCCAGGCCGCTCATTCTCTCCAGTACTCAAAAGCAAGCTTTGC;SVLEN=55;INSORD=0;SVID=sv0;SVTYPE=INS;SYMBOL=A 	GT	1|1
```

### DEL: A -> ''

POS gives the first deleted base.  END gives the last deleted base.  SVLEN gives the number of bases deleted.
The operation is marked in\_place in the IN\_PLACE field to indicate the source region was 
directly modified.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	6	sv0	N	<DEL>	100	PASS	END=60;OP_TYPE=NA;GRAMMAR=A->;VSET=0;SVLEN=55;SVID=sv0;SVTYPE=DEL;SYMBOL=A 	GT	1|1
```

### INV: A -> a

POS gives the first base of the inverted region.  END gives the last base of the inverted
region.  SVLEN gives the length of the inverted region.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	6	sv0	N	<INV>	100	PASS	END=60;OP_TYPE=NA;GRAMMAR=A->a;VSET=0;SVLEN=55;SVID=sv0;SVTYPE=INV;SYMBOL=A 	GT	1|1
```

### DUP: A -> AA

POS gives the first base of the duplicated region.  END gives the last base of the duplicated
region.  SVLEN gives the length of the duplicated region.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	6	sv0	N	<DUP>	100	PASS	END=60;OP_TYPE=NA;GRAMMAR=A->AA+;VSET=0;TARGET_CHROM=chr1;TARGET=61;SVLEN=55;SVID=sv0;SVTYPE=DUP;SYMBOL=A 	GT	1|1
```
## SVs represented as multiple VCF records

### rTRA: A_B -> B_A
A reciprocal translocation is represented by a combination of four records.
A is pasted at B's starting base, this is represented by the first record of type IDENTITY, the insertion target
is given by the TARGET and TARGET\_CHROM fields. A is also deleted from its location in the reference,
this operation is represented by the first record of type DEL.
The two other records represent the symmetrical operations affecting B.
All the records share the same GRAMMAR, SVID and SVTYPE.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1	sv0_0	N	<CUT-PASTE>	100	PASS	END=21;OP_TYPE=CUT-PASTE;GRAMMAR=A_B->B_A;VSET=0;TARGET_CHROM=chr1;TARGET=45;SVLEN=20;SVID=sv0;SVTYPE=rTRA;SYMBOL=A 	GT	1|1
chr1	45	sv0_1	N	<CUT-PASTE>	100	PASS	END=60;OP_TYPE=CUT-PASTE;GRAMMAR=A_B->B_A;VSET=0;TARGET_CHROM=chr1;TARGET=1;SVLEN=15;SVID=sv0;SVTYPE=rTRA;SYMBOL=B 	GT	1|1
```
## VCF fields
1. *CHROM*: Chromosome of the source region.
2. *POS*: Start position of the source region, 1 indexed.
3. *ID*: Unique identifier of the record containing the SV number. For SVs represented by multiple records, 
an additional suffix is appended to denote the specific record number within the group. 
4. *REF*: Reference allele for SNPs, N otherwise.
5. *ALT*: List of the alternative alleles for SNPs, SV type for SVs represented by a single record, operation type for SVs represented by multiple records.
6. *QUAL*: 100.
7. *FILTER*: PASS
8. *INFO*: Info fields described below.
9. *FORMAT*: GT
10. *SAMPLE*: Genotype 1|1 for homozygous variants, 1|0 or 0|1 for heterozygous ones.

### INFO fields
1. *END*: End of the source region, 1 indexed.
2. *OP_TYPE*: NA for SVs represented by a single record, the type of operation otherwise.
3. *GRAMMAR*: The corresponding SV grammar, not reported for tandem repeats.
4. *VSET*: Variant sets are numbered according to the order in which they appear in the config file.
5. *TARGET_CHROM*: Chromosome of the target point if the SV is not in place.
6. *TARGET*: Position of the target if the SV is not in place.
7. *SVLEN*: Length of the source region, corresponds to END - POS.
8. *INSORD*: Order in which the sequence corresponding to the record has been inserted. Used to disambiguate cases like A_ -> _AA, the first A on the right hand side will have INSORD=0 the second 1.
9. *SVID*: Unique identifier of the SV.
10. *SVTYPE*: Type of the SV.
11. *SYMBOL*: The letter from the grammar corresponding to the record source region. 

**Note:** For SNPs, VSET is the only reported INFO field.

## Output PAF file
The correct whole-genome alignment of the simulated sample haplotypes to the reference can optionally be
written out in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).
To enable this, add the setting `output_paf: True`.   The alignment will
be written to files named `sim.hapA.paf` and `sim.hapB.paf`, representing the alignment
to the reference of `sim.hapA.fa` and `sim.hapB.fa`, respectively.











