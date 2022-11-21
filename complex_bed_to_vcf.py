import argparse


def convert_to_vcf(bedfile, vcffile, template_vcf, include_simple):
    vcf_out = open(vcffile, 'w')

    with open(template_vcf, 'r') as templ:
        for line in templ:
            if line[0] == '#':
                vcf_out.write(line)
            else:
                break

    with open(bedfile, 'r') as bf:
        for line in bf:
            ln = line.split()
            chrm = ln[0]
            intervalA_start = int(ln[1])
            intervalA_end = int(ln[2])
            intervalB_start = int(ln[4])
            intervalB_end = int(ln[5])
            svtype = ln[-3]
            gt = ln[8]

            if svtype in ['dDUP', 'INV_dDUP']:
                vcf_out.write(f'{chrm}\t{intervalA_start}\t{svtype}_A\tN\t<{svtype}_A>\t.\tPASS\t'
                              f'END={intervalA_end};SVTYPE={svtype}_A;SVLEN={intervalA_end - intervalA_start}\tGT\t{gt}\n')
                vcf_out.write(f'{chrm}\t{intervalA_end}\t{svtype}_B\tN\t<{svtype}_B>\t.\tPASS\t'
                                 f'END={intervalB_end};SVTYPE={svtype}_B;SVLEN={intervalB_end - intervalA_end}\tGT\t{gt}\n')
            elif svtype == 'TRA':
                vcf_out.write(f'{chrm}\t{intervalA_start}\t{svtype}\tN\t<{svtype}>\t.\tPASS\t'
                              f'END={intervalA_end};SVTYPE={svtype};SVLEN={intervalA_end - intervalA_start}\tGT\t{gt}\n')
            elif svtype in include_simple:
                vcf_out.write(f'{chrm}\t{intervalA_start}\t{svtype}\tN\t<{svtype}>\t.\tPASS\t'
                              f'END={intervalA_end};SVTYPE={svtype};SVLEN={intervalA_end - intervalA_start}\tGT\t{gt}\n')

    vcf_out.close()    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert complex bedpe file to single-record vcf')
    parser.add_argument('--bed_file', type=str, help='Path to input bed file')
    parser.add_argument('--vcf_file', type=str, help='Path to output vcf file')
    parser.add_argument('--template_vcf', type=str, help='Template VCF to provide header info not given in bed file')
    parser.add_argument('--include_simple', nargs='*', help='List of simple types to allow in output vcf')
    args = parser.parse_args()

    convert_to_vcf(args.bed_file, args.vcf_file, args.template_vcf, args.include_simple)
