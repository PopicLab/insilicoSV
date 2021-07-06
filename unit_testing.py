from processing import FastaFile, Formater
import sys
import random

tests = {(10,15): "GAGTC",
         (67,70): "TAT",
         (140, 151): "AGTCCTATGAG",
         (145, 155): "Error"}

def test_fetch(tests):
    for (x,y) in tests:
        fasta = FastaFile("reference/test.txt")
        print(fasta)
        result = fasta.fetch(x,y)
        if result == tests[(x,y)]:
            print("Ok")
        else:
            print("Failed. Simulation result {} does not match true result {}".format(result, tests[(x,y)]))

def test_export():
    fasta = FastaFile("reference/test.txt")
    variants = {">Chromosome 19": [[3,3,"UUUUUUUUUU"]]}
    variants2 = {">Chromosome 21": [[4,10,"AAA"]]}
    fasta.reset_file("reference/test_out.txt")
    fasta.export_piece(variants, "reference/test_out.txt", verbose = True)
    fasta.export_piece(variants2, "reference/test_out.txt",  verbose = False)
    #fasta.close()

def test_target_pos():
    pass

output = ""
for x in range(10000):
    result = random.randint(1,4)
    if result == 1:
        output += "A"
    elif result == 2:
        output += "T"
    elif result == 3:
        output += "G"
    elif result == 4:
        output += "C"

print(output)