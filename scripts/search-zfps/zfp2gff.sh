#!/usr/bin/env bash

rg -f $1 \
    /Users/jonwells/Genomes/Mollusca/GCF_006345805.1_ASM634580v1_genomic.gff |
    rg -o -e "Parent=rna-(.+?);" |
    sed 's/Parent=rna-//' |
    sed 's/;//' |
    sort -u > parent_rnas.txt

rg -f parent_rnas.txt \
    /Users/jonwells/Genomes/Mollusca/GCF_006345805.1_ASM634580v1_genomic.gff \
    > ./octopus_transcripts.gff

rg '\tmRNA\t' ./octopus_transcripts.gff > ./octopus_transcripts_only.gff
rm parent_rnas.txt
