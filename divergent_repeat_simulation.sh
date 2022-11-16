# bash script to execute procedure of creating read set with divergent repeats
# inputs: - config files for R1 and R2 (the two references generated as part of the simulation procedure)

SCRIPT=$(realpath "$0")
SCRIPT_PATH=$(dirname "$SCRIPT")
SOURCE_REF=$1
CONFIG_R1=$2  # <-- insilicoSV randomized mode (optional SVs section input of - avoid_intervals: {vcf_path}
CONFIG_R2=$3  #<-- insilicoSV fixed mode config file that points to R2.vcf (or whatever it's called)
OUTPUT_PREFIX=$4  #<-- Typically just 'R'

# generate R1 and R2
python ${SCRIPT_PATH}/simulate.py ${SOURCE_REF} ${CONFIG_R1} ${OUTPUT_PREFIX}1
python ${SCRIPT_PATH}/simulate.py ${SOURCE_REF} ${CONFIG_R2} ${OUTPUT_PREFIX}2
# generate reads from hapA and hapB from the R2 fastas
dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H -R 0.30 -X 0.5 ${OUTPUT_PREFIX}2.hapA.fa ${OUTPUT_PREFIX}2.dwgsim.hap.0
dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H -R 0.30 -X 0.5 ${OUTPUT_PREFIX}2.hapB.fa ${OUTPUT_PREFIX}2.dwgsim.hap.1
# concat the R2 fastas
mv ${OUTPUT_PREFIX}2.dwgsim.hap.0.bfast.fastq.gz ${OUTPUT_PREFIX}2.dwgsim.hap.12.bfast.fastq.gz
cat ${OUTPUT_PREFIX}2.dwgsim.hap.1.bfast.fastq.gz >> ${OUTPUT_PREFIX}2.dwgsim.hap.12.bfast.fastq.gz
# map the R2 reads to R1 (either haplotype since we required div_dDUPs to be input to R1 as homozygous)
bwa index ${OUTPUT_PREFIX}1.hapA.fa
bwa mem -t20 -p ${OUTPUT_PREFIX}1.hapA.fa ${OUTPUT_PREFIX}2.dwgsim.hap.12.bfast.fastq.gz | samtools view -Sb - > ${OUTPUT_PREFIX}.bwamem.bam
samtools sort -@ 20 ${OUTPUT_PREFIX}.bwamem.bam -o ${OUTPUT_PREFIX}.bwamem.sorted.bam
samtools index -@ 20 ${OUTPUT_PREFIX}.bwamem.bam
# call position correcting script on R1 to make loci match with locations of repeats placed into genome
# --> we assume this script and locus adjustment script (fix_fiv_dDUP_vcf.py) are in same directory
python ${SCRIPT_PATH}/fix_div_dDUP_vcf.py --input_vcf ${OUTPUT_PREFIX}1.vcf --label_entire_event
# NB: flag --label_entire_event toggles whether the output vcf reports div. repeats in a single record, or
# uses two records to report interval A and interval a separately (including the flag sets behavior to the former)