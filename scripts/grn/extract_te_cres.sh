#!/usr/bin/env bash

###################################################################################################
## This script parses data from ENCODE iPS DF 19.11 cells and TRIM28 ChIP-exo experiments to
## identify candidate cis-regulatory elements (promoters and proximal enhancers) to identify
## regulatory relationships between KRAB-zinc finger proteins and transposable elements. As output
## produces three BED files:  
## 1. kzfp_candidate_promoters.bed
## 2. kzfp_candidate_proximal_enhancers.bed
## 3. kzfp_TRIM28_regions.bed
###################################################################################################

# Append associated KZFP gene labels to candidate CREs in intersected bedfile
function label_zfp_cre()
{
    infile=$1
    rg -o "gene=(.+?);" $infile | sed 's/gene=//' | sed 's/;//' > tmp_znf.txt
    paste -d '\t' $infile tmp_znf.txt |
    awk -F '\t' '{ 
            OFS="\t" 
        } {
        if ($13 ~ "PLS")
            print $10, $11, $12, $13 ";" $NF, $14, $7
        else
            print $10, $11, $12, $13 ";" $NF, $14, "."
        }'
    rm tmp_znf.txt
}

# Append associated TE labels to candidate CREs in intersected bedfile
function label_cre_te()
{
    infile=$1
    awk -F '\t' '{ 
            OFS="\t" 
        } {
            print $1, $2, $3, $4 ";" $10, $5, $6, $17
        }' $infile
}

# Append associated TE labels to TRIM28 peaks in intersected bedfile
function label_heterochromatin_te()
{
    infile=$1
    awk -F '\t' '{ 
            OFS="\t" 
        } {
            print $1, $2, $3, $4 ";" $9, $5, ".", $16
        }' $infile
}
###################################################################################################
## Get KZFP transcription start sites
###################################################################################################

awk -F '\t' '{
    OFS="\t"
    } {
    if ($7 == "+")
        print $1, $2, $3, $4, $4, $6, $7, $8, $9
    else
        print $1, $2, $3, $5, $5, $6, $7, $8, $9
    }' ../../data/gffs/kzfp_transcripts_only.gff |
    bedtools sort > ../../data/gffs/kzfp_transcripts_tss.gff

###################################################################################################
# Get KZFP promoters from iPS DF 19.11 cells using window of 500bp upstream and 300bp downstream of
# transcription start site. 
###################################################################################################

# Extract promoter-like sequences from cCRE file.
rg ';PLS' ../../data/cis-reg/ENCFF653NZY_T2T_ncbi.bed | bedtools sort > tmp_PLS.bed

# Define window within which to find promoters
bedtools slop \
    -i ../../data/gffs/kzfp_transcripts_tss.gff \
    -g ../../data/genome/GCF_009914755.1.genome \
    -l 500 \
    -r 300 \
    -s |
    bedtools sort > kzfp_promoter_region.bed

# Extract and label overlapping promoters 
bedtools intersect \
    -a kzfp_promoter_region.bed \
    -b tmp_PLS.bed \
    -wa \
    -wb > tmp.bed
label_zfp_cre tmp.bed > tmp2.bed

# label overlapping TEs
bedtools intersect \
    -a tmp2.bed \
    -b ~/Genomes/Mammalia/Hominoidea/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_rm.bed \
    -loj \
    -wa \
    -wb \
    -f 0.2 \
    -F 0.1 > tmp3.bed
label_cre_te tmp3.bed > ../../data/cis-reg/kzfp_candidate_promoters.bed

# Cleanup
rm tmp_PLS.bed tmp{,2,3}.bed kzfp_promoter_region.bed

###################################################################################################
# Get KZFP promoters from iPS DF 19.11 cells using window of 1kbp around TSS. 
###################################################################################################

# Extract proximal enhancer-like sequences from cCRE file.
rg ';pELS' ../../data/cis-reg/ENCFF653NZY_T2T_ncbi.bed | bedtools sort > tmp_pELS.bed

# Define window within which to find proximal enhancer-like regions
bedtools slop \
    -i ../../data/gffs/kzfp_transcripts_tss.gff \
    -g ../../data/genome/GCF_009914755.1.genome \
    -b 1000 |
    bedtools sort > kzfp_extended_promoter_region.bed

# Extract and label overlapping enhancers
bedtools intersect \
    -a kzfp_extended_promoter_region.bed \
    -b tmp_pELS.bed \
    -wa \
    -wb > tmp.bed
label_zfp_cre tmp.bed > tmp2.bed

# label overlapping TEs
bedtools intersect \
    -a tmp2.bed \
    -b ~/Genomes/Mammalia/Hominoidea/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_rm.bed \
    -loj \
    -wa \
    -wb \
    -f 0.2 \
    -F 0.1 > tmp3.bed
label_cre_te tmp3.bed  > ../../data/cis-reg/kzfp_candidate_proximal_enhancers.bed

# Cleanup
rm tmp_pELS.bed tmp{,2,3}.bed kzfp_extended_promoter_region.bed

###################################################################################################
# Get TRIM28-bound TEs within 10kb of KZFP TSS
###################################################################################################

# Define 10kb flanks either side of kzfp sequence.
bedtools slop \
    -i ../../data/gffs/kzfp_transcripts_only.gff \
    -g ../../data/genome/GCF_009914755.1.genome \
    -b 10000 |
    bedtools sort > kzfp_heterochromatin_region.bed

# Extract and label TRIM28-bound sites within 10kb of KZFPs
bedtools intersect \
    -a kzfp_heterochromatin_region.bed \
    -b ../../data/cis-reg/GSM2067350_KAP1_exo_H1_peaks_ncbi_T2T.bed \
    -wa \
    -wb > tmp.bed
label_zfp_cre tmp.bed |
    awk -F '\t' '{ OFS="\t" } { print $1, $2, $3, $4, $5 }' |
    bedtools sort > tmp2.bed

# label overlapping TEs
bedtools intersect \
    -a tmp2.bed \
    -b ~/Genomes/Mammalia/Hominoidea/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_rm.bed \
    -wa \
    -wb \
    -loj \
    -f 0.25 > tmp3.bed
label_heterochromatin_te tmp3.bed > ../../data/cis-reg/kzfp_TRIM28_regions.bed

# Cleanup
rm tmp{,2,3}.bed kzfp_heterochromatin_region.bed
