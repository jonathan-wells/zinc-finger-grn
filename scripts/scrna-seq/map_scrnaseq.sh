#!/usr/bin/env bash

# STAR version 2.7.10b
# Generate genome indices
/programs/STAR-2.7.10b/STAR \
    --genomeDir "../../data/expression/STAR-genome" \
    --genomeFastaFiles "/workdir/Genomes/Homo_sapiens/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna" \
    --runMode genomeGenerate \
    --runThreadN 12 \
    --sjdbGTFfile "/workdir/Genomes/Homo_sapiens/GCF_009914755.1/TEtranscripts.gtf" \
    --sjdbOverhang 74 \
    --limitSjdbInsertNsj 2000000
            

# Align reads
for readfile in $(cat ../../data/expression/SRR_Acc_List.txt); do
    echo Processing $readfile ...
    /programs/STAR-2.7.10b/STAR \
        --genomeDir "../../data/expression/STAR-genome" \
        --readFilesIn "../../data/expression/${readfile}.fastq" \
        --runThreadN 12 \
        --outSAMtype BAM Unsorted \
        --runMode alignReads \
        --outFilterMultimapNmax 250 \
        --winAnchorMultimapNmax 500 \
        --outMultimapperOrder Random \
        --alignIntronMax 500000 \
        --alignMatesGapMax 500000 \
        --outFileNamePrefix "../../data/expression/STAR-aligned/${readfile}_"
done
