#!/usr/bin/env bash

rg -f ../../data/hmmer-out/krab_zfp_names.txt \
    ~/Genomes/Mammalia/Hominoidea/GCF_009914755.1/genomic.gff |
    rg -o -e "Parent=rna-(.+?);" |
    sed 's/Parent=rna-//' |
    sed 's/;//' |
    sort -u > parent_rnas.txt

rg -f parent_rnas.txt \
    ~/Genomes/Mammalia/Hominoidea/GCF_009914755.1/genomic.gff \
    > ../../data/gffs/kzfp_transcripts.gff

rg '\tmRNA\t' ../../data/gffs/kzfp_transcripts.gff > ../../data/gffs/kzfp_transcripts_only.gff
rm parent_rnas.txt
