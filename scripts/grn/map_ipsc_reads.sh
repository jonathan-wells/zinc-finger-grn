#!/usr/bin/env bash

# STAR version 2.7.10b
# Generate genome indices
/programs/STAR-2.7.10b/STAR \
    --genomeDir "../../data/cis-reg/ipsc-expression/STAR-genome" \
    --genomeFastaFiles "/workdir/Genomes/Homo_sapiens/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna" \
    --runMode genomeGenerate \
    --runThreadN 32 \
    --sjdbGTFfile "/workdir/Genomes/Homo_sapiens/GCF_009914755.1/TEtranscripts.gtf" \
    --sjdbOverhang 74 \
    --limitSjdbInsertNsj 2000000
            

# Align reads
for readfile in $(cat ../../data/cis-reg/ipsc-expression/accession_list.txt); do
    echo Processing $readfile ...
    /programs/STAR-2.7.10b/STAR \
        --genomeDir "../../data/cis-reg/ipsc-expression/STAR-genome" \
        --readFilesIn "../../data/cis-reg/ipsc-expression/${readfile}_1.fastq.gz" "../../data/cis-reg/ipsc-expression/${readfile}_2.fastq.gz" \
        --readFilesCommand gunzip -c \
        --runThreadN 32 \
        --outSAMtype BAM Unsorted \
        --runMode alignReads \
        --outFilterMultimapNmax 250 \
        --winAnchorMultimapNmax 500 \
        --outMultimapperOrder Random \
        --alignIntronMax 500000 \
        --alignMatesGapMax 500000 \
        --outFileNamePrefix "../../data/cis-reg/ipsc-expression/STAR-aligned/${readfile}_"
done
