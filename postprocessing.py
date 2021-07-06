import os

def merge_fastq(file1, file2, output_file):
    '''
    file1 & file2: fastq file
    
    merges fastq files together so that ids will be unique
    '''
    with open(output_file, "w") as fout:
        with open(file1, "r") as f1:
            for line in f1:
                if line.startswith("@"):
                    fout.write(line.strip()[:-2] + "\n")
                else:
                    fout.write(line)

        with open(file2, "r") as f2:
            for line in f2:
                if line.startswith("@"):
                    fout.write(line.strip()[:-2] + "\n")
                else:
                    fout.write(line)

def prepare_for_bwa(reads1, reads2, output_files, root = None):
    if root:
        reads1 = [os.path.join(root, ele) for ele in reads1]
        reads2 = [os.path.join(root, ele) for ele in reads2]
        output_files = [os.path.join(root, ele) for ele in output_files]
    
    print(reads1, reads2)
    merge_fastq(reads1[0], reads1[1], output_files[0])
    merge_fastq(reads2[0], reads2[1], output_files[1])

root = "reads"
reads1 = ["test1.bwa.read1.fastq", "test2.bwa.read1.fastq"]
reads2 = ["test1.bwa.read2.fastq", "test2.bwa.read2.fastq"]
output_files = ["test.read1.fastq", "test.read2.fastq"]

prepare_for_bwa(reads1, reads2, output_files, root)

        