# Output representation

The details of the output representation for each built-in SV type are illustrated below.

Unless otherwise noted, the examples below are with respect to the following toy
reference:

```
>toyA
CCTAGTAGTAGCCCCTAGTAGTAGCCCCTAGTAGTAGCCCCTAGTAGTAGCCGAGAGAGATTT
```

## SVs represented as one VCF record

### INS: '' -> A

POS gives the reference base before which the novel sequence is inserted.
END is set equal to POS.   The inserted sequence is given in the INFO field INSSEQ.
The length of the inserted sequence is given as SVLEN.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
toyA	2	sv0	N	<INS>	100	PASS	END=2;SVTYPE=INS;VSET=0;INSSEQ=GTTCGATGTACTCCCAGGCCGCTCATTCTCTCCAGTACTCAAAAGCAAGCTTTGC;SVLEN=55;INSORD=0	GT	1|1
```

### DEL: A -> ''

POS gives the first deleted base.  END gives the last deleted base.  SVLEN gives the number of bases deleted.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
toyA	6	sv0	N	<DEL>	100	PASS	END=60;SVTYPE=DEL;VSET=0;SVLEN=55	GT	1|1
```

### INV: A -> A

POS gives the first base of the inverted region.  END gives the last base of the inverted
region.  SVLEN gives the length of the inverted region.

Example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
toyA	6	sv0	N	<INV>	100	PASS	END=60;SVTYPE=INV;VSET=0;SVLEN=55	GT	1|1
```

### DUP: A -> AA+

POS gives the first base of the duplicated region.  END gives the last base of the duplicated
region.  SVLEN gives the length of the duplicated region.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
toyA	6	sv0	N	<INV>	100	PASS	END=60;SVTYPE=INV;VSET=0;SVLEN=55	GT	1|1
```













