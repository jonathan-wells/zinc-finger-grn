#!/usr/bin/env bash

# Search for C2H2 ZF proteins
hmmsearch \
    --tblout '../../data/hmmer-out/zf_C2H2_PF00096.27.out' \
    '../../data/phmms/zf_C2H2_PF00096.27.hmm' \
    '/Users/jonwells/Genomes/Mammalia/Hominoidea/GCF_009914755.1/longest_protein_isoform.faa' &

hmmsearch \
    --tblout '../../data/hmmer-out/KRAB_PF01352.28.out' \
    '../../data/phmms/KRAB_PF01352.28.hmm' \
    '/Users/jonwells/Genomes/Mammalia/Hominoidea/GCF_009914755.1/longest_protein_isoform.faa' &

wait

rg -v '#' '../../data/hmmer-out/zf_C2H2_PF00096.27.out' |
    awk '{ print $1 }' |
    sort -u > zf_names.txt

rg -v '#' '../../data/hmmer-out/KRAB_PF01352.28.out' |
    awk '{ print $1 }' |
    sort -u > krab_names.txt

comm -12 zf_names.txt krab_names.txt > krab_zfp_names.txt
seqkit grep \
    -f krab_zfp_names.txt \
    '/Users/jonwells/Genomes/Mammalia/Hominoidea/GCF_009914755.1/longest_protein_isoform.faa' \
    > '../../data/seqs/krab_zfps.fa'

mv *_names.txt ../../data/hmmer-out/

