# Example SV Visualizations

Below are IGV pileup images of each of the predefined SV types with short read alignment and pair orientation highlighted 
to demonstrate the alignment signatures manifested by each SV.

### SNP

![SNP](gallery/SNP.png)

### INS
$\emptyset \rightarrow$ A

![INS](gallery/INS.png)

### DEL
A $\rightarrow \emptyset$

![DEL](gallery/DEL.png)

### DUP
A $\rightarrow$ AA+

![DUP](gallery/DUP.png)

### mCNV
A $\rightarrow$ A+

![mCNV](gallery/mCNV.png)

### INV\_DUP
A $\rightarrow$ Aa

![INV_DUP](gallery/INV_DUP.png)

### DUP\_INV
A $\rightarrow$ aa

![DUP_INV](gallery/DUP_INV.png)

### INV
A $\rightarrow$ a

![INV](gallery/INV.png)

### dDUP
A\_ $\rightarrow$ A\_A or \_A $\rightarrow$ A\_A

![dDUP](gallery/dDUP.png)

### INV\_dDUP
A\_ $\rightarrow$ A\_a or \_A $\rightarrow$ a\_A

![INV_dDUP](gallery/INV_dDUP.png)

### dDUP\_INV
A\_ $\rightarrow$ a\_a or \_A $\rightarrow$ a\_a

![dDUP_INV](gallery/dDUP_INV.png)

### nrTRA
A\_ $\rightarrow$ \_A or \_A $\rightarrow$ A\_

![nrTRA](gallery/nrTRA.png)

### rTRA
A\_B $\rightarrow$ B\_A or B\_A $\rightarrow$ A\_B

![rTRA](gallery/rTRA.png)

### dDUP_iDEL
A\_B $\rightarrow$ A\_A or B\_A $\rightarrow$ A\_A

![dDUP_iDEL](gallery/dDUP_iDEL.png)

### INS_iDEL
A\_B $\rightarrow$ \_A or B\_A $\rightarrow$ A\_

![INS\_iDEL](gallery/INS_iDEL.png)

### dupINV
AB $\rightarrow$ Aba

![dupINV](gallery/dupINV.png)

### INVdup
AB $\rightarrow$ baB

![INVdup](gallery/INVdup.png)

### delINV
AB $\rightarrow$ b

![delINV](gallery/delINV.png)

### delINVdel
ABC $\rightarrow$ b

![delINVdel](gallery/delINVdel.png)

### dupINVdup
ABC $\rightarrow$ AcbaC

![dupINVdup](gallery/dupINVdup.png)

### dupINVdel
ABC $\rightarrow$ Aba

![dupINVdel](gallery/dupINVdel.png)

### delINVdup
ABC $\rightarrow$ cbC

![delINVdup](gallery/delINVdup.png)

### Tandem Repeat Examples
### trCON
The motif is contracted, the number of repetitions is reduced by a value in 'repeat_count_change_range'.
![trCON](gallery/trCON.png)


### trEXP
The motif is expanded, 'repeat_count_change_range' repetitions are added.
![trEXP](gallery/trEXP.png)

### Context Aware Examples
Provides examples of the impact of the different placement options on the SV signatures.
### Contained Overlap
The DEL is contained in an Alu region.
![DEL_contained](gallery/DEL_contained.png)

### Containing Overlap
The DEL contains a L1 region.
![DEL_containing](gallery/DEL_containing.png)

### Partial Overlap
One and only one of the DEL breakends is in an L1HS region.
![DEL_partial](gallery/del_partial.png)

### Exact Overlap
The DEL breakends delimit an Alu region.
![DEL_exact](gallery/DEL_exact.png)