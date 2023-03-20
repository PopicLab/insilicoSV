# variation of divergent_repeat_simulation.sh to allow for adding divergent repeats onto a genome alongside
# other events

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
python ${SCRIPT_PATH}/div_repeat_postproc.py --div_dDUP_vcf ${OUTPUT_PREFIX}1.vcf --avoid_intervals
# --> output of the above: R1_avoid_intervals.vcf (to be included in the simple_edit config file)
# add the simple events to R2, avoiding the intervals taken up by the divergent repeats now inserted into R2
python ${SCRIPT_PATH}/simulate.py ${OUTPUT_PREFIX}2.hapA.fa ${CONFIG_SIMPLE_EDIT} ${OUTPUT_PREFIX}2_EDIT
# --> the vcf output of this will only specify the simple and dDUP events (we will use this in the merge step)
# generate reads from R2_EDIT (ref with div. repeats and simple events/dDUPs)
dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H -R 0.30 -X 0.5 ${OUTPUT_PREFIX}2_EDIT.hapA.fa ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.0
dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H -R 0.30 -X 0.5 ${OUTPUT_PREFIX}2_EDIT.hapB.fa ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.1
# concat the R2 fastqs
mv ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.0.bfast.fastq.gz ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.12.fastq.gz
cat ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.1.bfast.fastq.gz >> ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.12.fastq.gz
rm ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.1.bfast.fastq.gz
# map the R2 reads to R1 (either haplotype since we required div_dDUPs to be input to R1 as homozygous)
bwa index ${OUTPUT_PREFIX}1.hapA.fa
bwa mem -t20 -p ${OUTPUT_PREFIX}1.hapA.fa ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.12.fastq.gz | samtools view -Sb - > ${OUTPUT_PREFIX}.bwamem.bam
## -- remove fastq after bam is generated
rm ${OUTPUT_PREFIX}2_EDIT.dwgsim.hap.12.fastq.gz
samtools sort -@ 20 ${OUTPUT_PREFIX}.bwamem.bam -o ${OUTPUT_PREFIX}.bwamem.sorted.bam
## -- remove unsorted bam
rm ${OUTPUT_PREFIX}.bwamem.bam
samtools index -@ 20 ${OUTPUT_PREFIX}.bwamem.sorted.bam