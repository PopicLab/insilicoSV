# (recreating reference-editing script that I must have left on Hendrix)
SCRIPT=$(realpath "$0")
SCRIPT_PATH=$(dirname "$SCRIPT")
CONFIG=$1 #<- simulation config for simulate.py
FASTA_PATH=$2 #<- reference to be used as base for simulation
OUTPUT_PREFIX=$3 #<- filename prefix for output files (should include path to output dir)

python ${SCRIPT_PATH}/simulate.py ${FASTA_PATH} ${CONFIG} ${OUTPUT_PREFIX}
#/athena/ihlab/scratch/vpopic/software/LRSIM/dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H ${OUTPUT_PREFIX}.hapA.fa ${OUTPUT_PREFIX}.dwgsim.hap.0
#/athena/ihlab/scratch/vpopic/software/LRSIM/dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H ${OUTPUT_PREFIX}.hapB.fa ${OUTPUT_PREFIX}.dwgsim.hap.1
# --> should take dwgsim path as input in general case?
/athena/ihlab/scratch/crohlice/software/DWGSIM/dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H ${OUTPUT_PREFIX}.hapA.fa ${OUTPUT_PREFIX}.dwgsim.hap.0
/athena/ihlab/scratch/crohlice/software/DWGSIM/dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H ${OUTPUT_PREFIX}.hapB.fa ${OUTPUT_PREFIX}.dwgsim.hap.1
mv ${OUTPUT_PREFIX}.dwgsim.hap.0.12.fastq ${OUTPUT_PREFIX}.dwgsim.hap.12.fastq
cat ${OUTPUT_PREFIX}.dwgsim.hap.1.12.fastq >> ${OUTPUT_PREFIX}.dwgsim.hap.12.fastq
rm ${OUTPUT_PREFIX}.dwgsim.hap.1.12.fastq
bwa mem -t 20 -p ${FASTA_PATH} ${OUTPUT_PREFIX}.dwgsim.hap.12.fastq | samtools view -Sb - > ${OUTPUT_PREFIX}.bwamem.bam
rm ${OUTPUT_PREFIX}.dwgsim.hap.12.fastq
samtools sort -@ 20 ${OUTPUT_PREFIX}.bwamem.bam > ${OUTPUT_PREFIX}.bwamem.sorted.bam
rm ${OUTPUT_PREFIX}.bwamem.bam
samtools index -@ 20 ${OUTPUT_PREFIX}.bwamem.sorted.bam
