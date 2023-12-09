#!/usr/bin/env python3

import sys
from collections import defaultdict


def parse_cisbp_metdata():
    """Reads CisBP metadata and returns most recent motif for each KZFP"""
    cisbp_dir = '../../data/cis-reg/CisBP_Homo_sapiens_2023_12_08_10_56_am'
    # Extract list if KZFPs
    with open('../../data/cis-reg/kzfp_list.txt') as input:
        kzfp_list = [line.strip() for line in input] 

    # Parse CisBP metadata file
    motifdict = defaultdict(list)
    with open(f'{cisbp_dir}/TF_Information.txt') as input:
        header = input.readline().strip('\n').split('\t')
        for line in input:
            line = line.strip('\n').split('\t')
            motif, kzfp, source, year = line[3], line[6], line[16], line[18]
            if motif == '.':
                continue
            year = int(year[:4])
            if kzfp not in kzfp_list:
                continue
            motifdict[kzfp].append((motif, source, year))

    # Select most recent motif and write to stdout
    for kzfp, val in motifdict.items():
        motif = max(val, key=lambda x: x[2])[0]
        sys.stdout.write(f'{kzfp}\t{motif}\n')

if __name__ == '__main__':
    parse_cisbp_metdata()
