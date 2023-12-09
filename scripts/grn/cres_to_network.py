#!/usr/bin/env python3

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
            else:
                if te != '.':
                    znf_to_te[znf].add(te)
    return znf_to_te

def write_edges(znf_to_te, outfile):
    """Convert znf_to_te dictionary into tab-separated list of KZFP-TE pairs."""
    with open(outfile, 'w') as output:
        for kzfp, val in znf_to_te.items():
            if len(val) == 0:
                output.write(f'{kzfp}\n')
            else:
                for te in val:
                    output.write(f'{kzfp}\t{te}\n')



def main():
    """Write paired edges for promoters, proximal enhancers and TRIM28 peaks."""
    trim28_te = parse_cisreg_bed('../../data/cis-reg/cCRE-bed/kzfp_TRIM28_regions.bed', trim28=True)
    promoter_te = parse_cisreg_bed('../../data/cis-reg/cCRE-bed/kzfp_candidate_promoters.bed', trim28=False)
    prox_enhancer_te = parse_cisreg_bed('../../data/cis-reg/cCRE-bed/kzfp_candidate_proximal_enhancers.bed', trim28=False)

    write_edges(trim28_te, '../../data/cis-reg/kzfp_TRIM28_regions.txt')
    write_edges(promoter_te, '../../data/cis-reg/kzfp_candidate_promoters.txt')
    write_edges(prox_enhancer_te, '../../data/cis-reg/kzfp_candidate_proximal_enhancers.txt')

if __name__ == '__main__':
    main()
