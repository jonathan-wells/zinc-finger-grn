#!/usr/bin/env python3

import re

def parse_cisreg_bed(bedfile, trim28=False):
    """Read KZFP-CRE BED file and return as dictionary from KZFPs to assoc. TEs.
    
    Arguments:
        bedfile - BED format file generated from `extract_te_cres.sh` script
        trim28 - Boolean value describing whether BED file contains TRIM28 peaks.
    Returns:
        znf_to_te - dictionary mapping KZFPs to set of associated TEs.
    """
    znf_to_te = {}
    with open(bedfile) as infile:
        for line in infile:
            line = line.strip('\n').split('\t')
            
            if trim28:
                trim28_peak, znf, te = line[3].split(';')
            else:
                encode_id, cre_type, znf, te = line[3].split(';')
            if znf not in znf_to_te.keys():
                znf_to_te[znf] = set()
            if te != '.':
                znf_to_te[znf].add(te)
    return znf_to_te

def parse_kzfp_targets(qthresh):
    with open('../../data/cis-reg/enrich_kzfp_perSubfam/zfp_file_list.txt') as input:
        kzfps = {line.split('\t')[0]: line.strip().split('\t')[1] for line in input}
    
    kzfp_to_te = {}
    for kzfp, filename in kzfps.items():
        kzfp_to_te[kzfp] = set()
        with open(f'../../data/cis-reg/enrich_kzfp_perSubfam/{filename}') as input:
            header = input.readline()
            for line in input:
                line = line.strip().split('\t')
                te, qvals = line[0], [float(x) for x in line[1:4]]
                if te == 'nonTE':
                    break
                # if qvals[0] > qthresh and qvals[1] > qthresh and qvals[2] > qthresh:
                if qvals[2] > qthresh:
                    continue
                kzfp_to_te[kzfp].add(te)
    return kzfp_to_te

def parse_tf_te_fimo(fimofile, qthresh):
    tf_to_target = {}
    with open(fimofile) as input:
        header = input.readline()
        for line in input:
            if line == '\n':
                break
            line = line.strip().split('\t')
            tf = line[1]
            target = line[2]
            qvalue = float(line[8])
            if qvalue > qthresh:
                continue
            if tf not in tf_to_target:
                tf_to_target[tf] = set()
            tf_to_target[tf].add(target)
        return tf_to_target

def parse_tf_kzfp_fimo(fimofile, qthresh):
    tf_to_target = {}
    zfp_pattern = re.compile(r'[\w-]+;PLS;(.+)?;')
    with open(fimofile) as input:
        header = input.readline()
        for line in input:
            if line == '\n':
                break
            line = line.strip().split('\t')
            tf = line[1]
            target = line[2]
            qvalue = float(line[8])
            if qvalue > qthresh:
                continue
            kzfphit = re.match(zfp_pattern, target)
            if kzfphit:
                target = kzfphit.group(1)
            if tf not in tf_to_target:
                tf_to_target[tf] = set()
            tf_to_target[tf].add(target)
        return tf_to_target

def write_edges(znf_to_te, outfile):
    """Convert znf_to_te dictionary into tab-separated list of KZFP-TE pairs."""
    with open(outfile, 'w') as output:
        for kzfp, val in znf_to_te.items():
            if len(val) == 0:
                output.write(f'{kzfp}\n')
            else:
                for te in val:
                    output.write(f'{kzfp}\t{te}\n')


def build_final_edge_list():
    """Takes output from various sources to compile list of edges for ZF Network"""
    tf2te = parse_tf_te_fimo('../../data/cis-reg/te-fimo/fimo.tsv', 0.01)
    tf2kzfp = parse_tf_kzfp_fimo('../../data/cis-reg/kzfp-fimo/fimo.tsv', 0.01)
    kzfp2te = parse_kzfp_targets(0.0001)
    te2kzfp = {}
    for kzfp, tes in kzfp2te.items():
        for te in tes:
            if te not in te2kzfp:
                te2kzfp[te] = set([kzfp])
            else:
                te2kzfp[te].add(kzfp)
    
    trim28_te = parse_cisreg_bed('../../data/cis-reg/cCRE-bed/kzfp_TRIM28_regions.bed', trim28=True)
    kzfp2kzfp = {}
    for kzfp_b, tes in trim28_te.items():
        for te in tes:
            if te not in te2kzfp:
                continue
            kzfps_a = te2kzfp[te]
            for kzfp_a in kzfps_a:
                if kzfp_a not in kzfp2kzfp:
                    kzfp2kzfp[kzfp_a] = set([kzfp_b])
                else:
                    kzfp2kzfp[kzfp_a].add(kzfp_b)
    with open('../../data/cis-reg/final_edge_list.txt', 'w') as output:
        for dictionary in (tf2te, tf2kzfp, kzfp2te, kzfp2kzfp):
            for key, values in dictionary.items():
                for val in values:
                    output.write(f'{key}\t{val}\n')
     
def main():
    build_final_edge_list()

if __name__ == '__main__':
    main()
