#!/usr/bin/env bash

# STAR version 2.7.11a
# Generate genome indices
STAR \
    --genomeDir "../../data/cis-reg/ipsc-expression/STAR-genome" \
    --genomeFastaFiles "/Users/jonwells/Genomes/Mammalia/Hominoidea/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna" \
    --runMode genomeGenerate \
    --runThreadN 8 \
    --sjdbGTFfile "/Users/jonwells/Genomes/Mammalia/Hominoidea/GCF_009914755.1/TEtranscripts.gtf" \
    --sjdbOverhang 74

# Align reads
for readfile in $(cat ../../data/cis-reg/ipsc-expression/accession_list.txt); do
    echo Processing $readfile ...
    STAR \
        --genomeDir "../../data/cis-reg/ipsc-expression/STAR-genome" \
        --readFilesIn "../../data/cis-reg/ipsc-expression/${readfile}_1.fastq.gz" "../../data/cis-reg/ipsc-expression/${readfile}_2.fastq.gz" \
        --readFilesCommand gunzip -c \
        --runThreadN 8 \
        --outSAMtype BAM Unsorted \
        --runMode alignReads \
        --outFilterMultimapNmax 250 \
        --winAnchorMultimapNmax 500 \
        --outMultimapperOrder Random \
        --alignIntronMax 500000 \
        --alignMatesGapMax 500000 \
        --outFileNamePrefix "../../data/cis-reg/ipsc-expression/STAR-aligned/${readfile}_"
done
