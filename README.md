# insilicoSV
insilicoSV is a software to design and simulate complex structural variants, both novel and known. 

## Requirements  (Prerequisites)
* Python 3.5 and up [Install](https://link-for-setup-guide)

## Installation

`$ pip install -r requirements.txt`

## To Run
`python simulate.py <ref.fna> <par.yaml> <hap1.fna> <hap2.fna> <out.bed>`

## Usage

insilicoSV takes in two input files: the reference genome and a yaml configuration file. Following the simulation, it outputs two haplotype files and a BEDPE file. 

### Reference Genome
This file should be in FASTA format.

### Parameter File
The configuration yaml file specifies the range of the lengths of the SVs along with the number to simulate. All configurations for structural variants should be put under the "SVs" key. For each SV, the following parameters are available:
1. *type*: str, insilicoSV supports a predefined list of SVs and allows users to enter a custom transformation. Either "Custom" or one of the 16 predefined SV types named in the below table should be entered.
2. *number*: int, describes how many of the specified SV type to simulate
3. *min_length*: int or list, if an integer is provided, insilcoSV assumes that each event's length within a SV must fall between the min_length and max_length. Entering a list offers customization by specifying a different range for each event.
4. *max_length*: int or list, must be the same type as min_length, note that max_length >= min_length >= 0 for all elements in each
5. *source [optional]*: Source sequence for a custom SV, see below to find instructions on how to create a transformation
6. *target [optional]*: Target sequence for a custom SV, see below to find instructions on how to create a transformation

Please see the table and picture for the list of predefined classes.

![Alt text](imgs/complex_sv_classes_diagram.webp)

| SV Type | Transformation |
|---------|----------------|
| INS | "" -> "A" |
| DEL | "A" -> "" |
| INV | "A" -> "a" |
| DUP | "A" -> "AA'" |
| TRA | "A_B" -> "B_A" |
| dupINVdup | "ABC" -> "Ac'ba'C" |
| delINVdel | "ABC" -> "b" |
| delINVdup | "ABC" -> "c'bC" |
| dupINVdel | "ABC" -> "Aba'" |
| delINV | "AB" -> "b" |
| INVdel | "AB" -> "a" |
| dDUP-iDEL | "A_B" -> "A_A'" |
| INS-iDEL | "A_B" -> "_A" |
| dupINV | "AB" -> "Aba'" |
| INVdup | "AB" -> "b'aB" |
| dDUP | "A_" -> "A_A'" |

A custom SV consists of a user-generated transformation with a source and target sequence of "symbols," most of which are alphabetical letters. Some examples of a source would be "ABC" and "A_B_C," while some examples of the target would be "a'AB" or "A_b'Bc'_C." 

insilicoSV maps a randomly generated piece of the reference to each of the symbols in the source and recompiles the affected region with the target. For instance, "AB -> "A" would remove the fragment marked as symbol B. 

Duplications of the original symbol must be marked with a "'", while inversions are marked in lowercase. Therefore, a transformation "ABC" -> "ABA'c" would duplicate "A" after "B" and invert the fragment C. 

All symbols in the source sequence MUST be unique to create a one-to-one mapping between symbol and reference fragment.

### BEDPE File

Each line/entry will have the following parameters:
1. *source_chr*: The source chromosome of the event
2. *source_start*: Start position on the source_chr, zero-based indexing
3. *source_end*: End position on the source_chr, one-based indexing
4. *target_chr*: The target chromosome of the event
5. *target_start*: Start position on the target chr, zero-based indexing
6. *target_end*: End position on the target chr, one-based indexing
7. *event_type*: Describes the transformation made by the event, either an INS, DEL, INV, TRA, DUP, INVDUP, or INVTRA. Dispersed duplications have an attached "d" at the front.
8. *event_size*: Size of the fragment impacted by the event on the source chromosome 
9. *parent_sv*: Describes the parent SV the event is a component of, for instance "dupINVdup." If a custom SV was provided, the name becomes "source>target"
10. *nth_sv*: int, index to count up each SV (note: not the events). All events of a SV belong in the same index.
11. *order*: int, for insertion-like operations such as TRA, INS, or DUP, the "order" index describes in which order events that target the same position were compiled. Events with INV and DEL have an order of 0.

## Usage example


## How to Contribute
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change. Please make sure to update tests as appropriate. If you'd like to contribute, please fork the repository and make changes as you'd like. Pull requests are warmly welcome.

Steps to contribute:
1. Fork this repository (link to your repository)
2. Create your feature branch (git checkout -b feature/insilicoSV)
3. Commit your changes (git commit -am 'Add some feature')
4. Push to the branch (git push origin feature/insilicoSV)
5. Create a new Pull Request
