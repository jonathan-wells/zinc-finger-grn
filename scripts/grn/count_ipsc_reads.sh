#!/usr/bin/env bash

while read sample; do
    echo $sample
    /programs/TEtranscripts-2.2.3/bin/TEcount \
        -b "../../data/cis-reg/ipsc-expression/STAR-aligned/${sample}_Aligned.out.bam" \
        --GTF /workdir/Genomes/Homo_sapiens/GCF_009914755.1/genomic.gtf \
        --TE /workdir/Genomes/Homo_sapiens/GCF_009914755.1/TEtranscripts_rmsk.gtf \
        --mode multi \
        --format BAM \
        --stranded reverse \
        --project "../../data/cis-reg/ipsc-expression/TEcount-out/${sample}" 
done < ../../data/cis-reg/ipsc-expression/accession_list.txt

