# variation of divergent_repeat_simulation.sh to allow for adding divergent repeats onto a genome alongside
# simple events and dispersed dups

SCRIPT=$(realpath "$0")
SCRIPT_PATH=$(dirname "$SCRIPT")
SOURCE_REF=$1
CONFIG_R1=$2
CONFIG_R2=$3
OUTPUT_PREFIX=$4
CONFIG_SIMPLE_EDIT=$5  #<-- (...and input divergent repeat intervals in SVs section input with - avoid_intervals: {vcf_path})

# generate R1 and R2
python ${SCRIPT_PATH}/simulate.py ${SOURCE_REF} ${CONFIG_R1} ${OUTPUT_PREFIX}1
python ${SCRIPT_PATH}/simulate.py ${SOURCE_REF} ${CONFIG_R2} ${OUTPUT_PREFIX}2
# create the vcf whose intervals represent A and A' for all events in R2
python ${SCRIPT_PATH}/fix_div_dDUP_vcf.py --div_dDUP_vcf ${OUTPUT_PREFIX}1.vcf --avoid_intervals
# --> output of the above: R1_avoid_intervals.vcf (to be included in the simple_edit config file)
# add the simple events to R2, avoiding the intervals taken up by the divergent repeats now inserted into R2
python ${SCRIPT_PATH}/simulate.py ${OUTPUT_PREFIX}2.hapA.fa ${CONFIG_SIMPLE_EDIT} ${OUTPUT_PREFIX}2_EDIT
# --> NB: the vcf output of this will only specify the simple and dDUP events (we will use this in the merge step)
# generate reads from R2_EDIT (ref with div. repeats and simple events/dDUPs)
#/athena/ihlab/scratch/vpopic/software/LRSIM/dwgsim
/athena/ihlab/scratch/crohlice/software/DWGSIM/dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H -R 0.30 -X 0.5 ${OUTPUT_PREFIX}2_EDIT.hapA.fa ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.0
/athena/ihlab/scratch/crohlice/software/DWGSIM/dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H -R 0.30 -X 0.5 ${OUTPUT_PREFIX}2_EDIT.hapB.fa ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.1
# concat the R2 fastqs
mv ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.0.12.fastq ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.12.fastq
cat ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.1.12.fastq >> ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.12.fastq
rm ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.1.12.fastq
# map the R2 reads to R1 (either haplotype since we required div_dDUPs to be input to R1 as homozygous)
bwa index ${OUTPUT_PREFIX}1.hapA.fa
bwa mem -t20 -p ${OUTPUT_PREFIX}1.hapA.fa ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.12.fastq | samtools view -Sb - > ${OUTPUT_PREFIX}.bwamem.bam
## -- remove fastq after bam is generated
rm ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.12.fastq
samtools sort -@ 20 ${OUTPUT_PREFIX}.bwamem.bam -o ${OUTPUT_PREFIX}.bwamem.sorted.bam
## -- remove unsorted bam
rm ${OUTPUT_PREFIX}.bwamem.bam
samtools index -@ 20 ${OUTPUT_PREFIX}.bwamem.sorted.bam
## optional step to create merged vcf with decoy div. repeats and non-decoy events listed together
#python ${SCRIPT_PATH}/fix_div_dDUP_vcf.py --div_dDUP_vcf ${OUTPUT_PREFIX}1.vcf --merge_input_vcf ${OUTPUT_PREFIX}2_EDIT_SIMPLE_DOUBLELABEL.vcf --merge_output_vcf ${OUTPUT_PREFIX}1_MERGED.vcf --label_entire_event
## NB: flag --label_entire_event toggles whether the output vcf reports div. repeats in a single record, or
## uses two records to report interval A and interval a separately (including the flag sets behavior to the former)