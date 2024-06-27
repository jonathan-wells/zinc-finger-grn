#!/usr/bin/env bash

c=1
while read sample; do
    echo $sample
    /programs/TEtranscripts-2.2.3/bin/TEcount \
        -b "../../data/expression/STAR-aligned/${sample}_Aligned.out.bam" \
        --GTF /workdir/Genomes/Homo_sapiens/GCF_009914755.1/genomic.gtf \
        --TE /workdir/Genomes/Homo_sapiens/GCF_009914755.1/TEtranscripts_rmsk.gtf \
        --mode multi \
        --format BAM \
        --stranded reverse \
        --project "../../data/expression/TEcount-out/${sample}" & 
    if [ $((c % 10)) -eq 0 ]; then
        wait
    fi
    c=$((c + 1))
done < ../../data/expression/SRR_Acc_List.txt

