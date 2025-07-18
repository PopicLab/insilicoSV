import argparse
from collections import defaultdict

def complement_bed(bed, complement_path=None):
    if not complement_path:
        complement_path = bed.rsplit('.', 1)[0] + '_complement.bed'

    regions = defaultdict(list)
    with open(bed, 'r') as f:
        for line in f.readlines():
            line = line.strip().split()
            regions[line[0]].append(line)

    with open(complement_path, 'w') as f:
        for chrom, chrom_regions in regions.items():
            chrom_regions.sort(key=lambda x: int(x[1]))
            for idx, region in enumerate(chrom_regions[1:]):
                previous_region = chrom_regions[idx]
                if int(previous_region[2]) > int(region[1]):
                    raise(Exception(f'Region {region} overlaps previous region {previous_region}. '
                                    f'Complement is not defined for overlapping regions.'))
                if previous_region[2] == region[1]: continue
                f.write('\t'.join([region[0], previous_region[2], region[1], 'complement']) + '\n')

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-b', '--bed', required=True, help='Input BED file to complement')
    parser.add_argument('-o', '--output', required=False, help='Complemented output BED file path')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    complement_bed(args.bed, args.output)