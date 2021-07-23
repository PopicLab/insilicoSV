import argparse
import os

def collect_args():
    parser = argparse.ArgumentParser(description='insilicoSV is a software to design and simulate complex structural variants, both novel and known.')
    parser.add_argument("ref", help="FASTA reference file")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("hap1", help="First output FASTA file (first haplotype)")
    parser.add_argument("hap2", help = "Second output FASTA file (second haplotype)")
    parser.add_argument("bedpe", help = "BEDPE file location to store variant info")
    parser.add_argument("-r", "--root", action="store", metavar="DIR", dest="root_dir", help="root directory for all files given as positional arguments")

    args = parser.parse_args()
    io = [args.ref, args.config, args.hap1, args.hap2, args.bedpe]
    if args.root_dir:
        io = [os.path.join(args.root_dir, ele) for ele in io]
    return io