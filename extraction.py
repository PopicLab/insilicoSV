from os import truncate


with open("reference/ref_genome.fna","r") as fin:
    with open("reference/ref_chr21.fna","w") as fout:
        go = False
        for line in fin:
            if line.startswith(">NC_000021.9"):
                print("Chromosome 21 Detected")
                go = True
            elif line.startswith(">") and go:
                go = False
                break
        
            if go:
                fout.write(line)


id = ">NC_000021.9 Homo sapiens chromosome 21, GRCh38.p13 Primary Assembly"

'''with open("reference/ref_genome.fna","r") as fin:
    for line in fin:
        if line.startswith(">"):
            print(line)
'''