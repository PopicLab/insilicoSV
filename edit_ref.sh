# (recreating reference-editing script that I must have left on Hendrix)
SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
CONFIG=$1 #<- simulation config for simulate.py
FASTA_PATH=$2 #<- reference to be used as base for simulation
OUTPUT_PREFIX=$3 #<- filename prefix for output files (should include path to output dir)

python ${SCRIPTPATH}/simulate.py ${FASTAPATH} ${CONFIG} ${OUTPUT_PREFIX}
dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H ${OUTPUT_PREFIX}.hapA.fa ${OUTPUT_PREFIX}.dwgsim.hap.0
dwgsim -C 30 -1 151 -2 151 -y 0 -S 0 -c 0 -m /dev/null -H ${OUTPUT_PREFIX}.hapB.fa ${OUTPUT_PREFIX}.dwgsim.hap.1
mv ${OUTPUT_PREFIX}.dwgsim.hap.0.bfast.fastq.gz ${OUTPUT_PREFIX}.dwgsim.hap.12.bfast.fastq.gz
cat ${OUTPUT_PREFIX}.dwgsim.hap.1.bfast.fastq.gz >> ${OUTPUT_PREFIX}.dwgsim.hap.12.bfast.fastq.gz
rm ${OUTPUT_PREFIX}.dwgsim.hap.1.bfast.fastq.gz
bwa mem -t20 -p ${FASTA_PATH} ${OUTPUT_PREFIX}.dwgsim.hap.12.bfast.fastq.gz | samtools view -Sb - > ${OUTPUT_PREFIX}.bwamem.bam
rm ${OUTPUT_PREFIX}.dwgsim.hap.12.bfast.fastq.gz
samtools sort -@ 20 ${OUTPUT_PREFIX}.bwamem.bam > ${OUTPUT_PREFIX}.bwamem.sorted.bam
rm ${OUTPUT_PREFIX}.bwamem.bam
samtools index -@ 20 ${OUTPUT_PREFIX}.bwamem.sorted.bam
