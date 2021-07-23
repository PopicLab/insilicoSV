# insilicoSV
insilicoSV is a software to design and simulate complex structural variants, both novel and known. 

## To Run
`python simulate.py <ref.fna> <par.yaml> <hap1.fna> <hap2.fna> <out.bed>`

## Parameter File

The configuration yaml file specifies the range of the lengths of the SVs along with the number to simulate. All SVs are distributed in a roughly even manner. Each element must have four attributes: type of SV (use below table for reference), number, minimum length, and maximum length, where max_length >= min_length >= 0. The values for min_length and max_length may either be an array or integer. If it is an integer, insilicoSV chooses a random number within the given range to reflect the size of *each* event/letter in the SV. If it is an array, the corresponding integers between the lists represent the size range for the event of the same position (e.g the first integers of both arrays give the range for the first event of the SV, and so on). Note that if the value is a list, the number of elements must match the number of events in the chosen SV. 

![Alt text](imgs/complex_sv_classes_diagram.webp)

For the configuration file, use this key when labelling the structural variant:

1 = Insertion\
2 = Deletion\
3 = Inversion\
4 = Duplication\
5 = Translocation\
6 = dupINVdup\
7 = delINVdel\
8 = delINVdup\
9 = dupINVdel\
10 = delINV\
11 = INVdel\
12 = dDUP-iDEL\
13 = INS-iDEL\
14 = dupINV\
15 = INVdup\
16 = dDUP
