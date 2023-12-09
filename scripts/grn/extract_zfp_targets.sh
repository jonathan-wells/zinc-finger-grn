#!/usr/bin/env bash

####################################################################################################
## This script KZFP-KZFP interactions by directly scanning for KZFP binding motifs within TRIM28 
## ChIP-exo peaks. It does so for both peaks associated with a particular TE, and for peaks without
## associated TEs.
####################################################################################################

# Extract fasta file of KZFP-associated TRIM28 peaks
bedtools getfasta \
    -fi ~/Genomes/Mammalia/Hominoidea/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
    -bed ../../data/cis-reg/cCRE-bed/kzfp_TRIM28_regions.bed \
    -name > ../../data/cis-reg/kzfp_TRIM28_regions.fa

# Extract background letter freqs for meme
/usr/local/meme-5.3.3/meme/libexec/meme-5.3.3/fasta-get-markov \
    -dna \
    ~/Genomes/Mammalia/Hominoidea/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
    genome_background.bg

# Clear pre-exisiting motif file
if [ -f ../../data/cis-reg/kzfp_motifs.meme ]; then
    rm ../../data/cis-reg/kzfp_motifs.meme
fi

# Convert CisBP PWMs to meme format motif files.
CISBPDIR="../../data/cis-reg/CisBP_Homo_sapiens_2023_12_08_10_56_am"
while read -r line; do
    echo $kzfp
    kzfp=$(awk '{ print $1 }' <<< $line) 
    motif=$(awk -v path="$CISBPDIR" '{ print path "/pwms_all_motifs/" $2 ".txt" }' <<< $line)
    awk 'NR > 1 { print $2, $3, $4, $5 }' $motif > tmp.pwm
    /usr/local/meme-5.3.3/scripts/matrix2meme -bg genome_background.bg < tmp.pwm |
        sed "s/MOTIF 1/MOTIF $kzfp/" >> ../../data/cis-reg/kzfp_motifs.meme
done < <(./parse_cisbp.py)

# # Run FIMO to identify KZFPs targeting KZFP-TRIM28 peaks
fimo \
    -oc ../../data/cis-reg/kzfp-trim28-fimo \
    -bfile genome_background.bg \
    ../../data/cis-reg/kzfp_motifs.meme \
    ../../data/cis-reg/kzfp_TRIM28_regions.fa

fimo \
    -oc ../../data/cis-reg/kzfp-trim28-noTE-fimo \
    -bfile genome_background.bg \
    ../../data/cis-reg/kzfp_motifs.meme \
    ../../data/cis-reg/kzfp_TRIM28_regions_noTE.fa

rm tmp.pwm
rm genome_background.bg
