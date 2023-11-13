#!/usr/bin/env bash

# VERY UNFINISHED

# TODO: Extend this to search in multiple species.

# Search for C2H2 ZF proteins
hmmsearch \
    --tblout '../../data/hmmer-out/zf_C2H2_PF00096.27.out' \
    '../../data/phmms/zf_C2H2_PF00096.27.hmm' \
    '/Users/jonwells/Genomes/Mollusca/GCF_006345805.1_ASM634580v1_longest_isoform.faa' 

# If number of ZNF hits matching reporting threshold (col 17) is reached, extract label
rg -v '#' '../../data/hmmer-out/zf_C2H2_PF00096.27.out' |
    awk '{ if ( $17 >= 5 ) print $1 }' |
    sort -u > zf_names.txt

seqkit grep \
    -f zf_names.txt \
    '/Users/jonwells/Genomes/Mollusca/GCF_006345805.1_ASM634580v1_longest_isoform.faa' \
    > '../../data/seqs/octopus_sinensis_zfps.fa'


