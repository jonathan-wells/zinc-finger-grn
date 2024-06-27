#!/usr/bin/env bash

####################################################################################################
## This script uses experimentally determined DNA-binding motifs from the JASPAR 2022 core
## non-redundant database to identify transcription factors potentially associated with KZFPs and
## TEs in the human genome. To reduce the search space, we limit the TFs used to those with at least
## 10 reads in an iPSC dataset. Similarly, we limit the TEs searched to those which are
## significantly targeted by at least one KZFP.
####################################################################################################

MINCOUNTS=10  # Minimum number of reads to include TFs in eventual network
PADTHRESH=1e-5  # adjusted p-value threshold for including TE in network 

# Extract list of TFs represented in JASPAR core
rg "MOTIF" ../../data/cis-reg/JASPAR2022_CORE_non-redundant_pfms_meme/*.meme |
    awk '{ OFS="\t" } { print $2, $3 }'> motif_dict.txt
awk '{ print "gene-"$2}' motif_dict.txt | sed 's/$/"/' | sed 's/^/"/' > search_names.txt

# Extract TFs with min number of reads
rg --ignore-case -f search_names.txt ../../data/cis-reg/ipsc-expression/TEcount-out/ERR947017.cntTable |
    awk -v mc=$MINCOUNTS '{ if ($2 >= mc) print $1 }' |
    sed 's/"/\\s/' | sed 's/gene-//' | sed 's/"/$/' > mincount_names.txt 

# Filter out KZFPs from list of TF-motifs to scan.
sed 's/^/\\s/' ../../data/cis-reg/kzfp_list.txt | 
    sed 's/$/$/' ../../data/cis-reg/kzfp_list.txt > omit_kzfps.txt
rg --ignore-case -f mincount_names.txt motif_dict.txt |
    rg -v --ignore-case -f omit_kzfps.txt |
    awk '{ print "../../data/cis-reg/JASPAR2022_CORE_non-redundant_pfms_meme/"$1".meme" }' \
    > meme_names.txt
cat $(cat meme_names.txt) > ../../data/cis-reg/ipsc_motifs.meme

# Get genome background for motif scanning
/usr/local/meme-5.3.3/meme/libexec/meme-5.3.3/fasta-get-markov \
    -dna \
    ~/Genomes/Mammalia/Hominoidea/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
    genome_background.bg

# Extract TEs targeted by at least one KZFP and combine with KZFP candidate promoters
if [ -f te_fimo_sequences.txt ]; then
    rm te_fimo_sequences.txt
fi
if [ -f kzfp_fimo_sequences.txt ]; then
    rm kzfp_fimo_sequences.txt
fi

for file in $(ls ../../data/cis-reg/enrich_kzfp_perSubfam/*te_subfam.txt); do
    awk -v p="$PADTHRESH" '{ if (($2 <= p) && ($3 < p) && ($4 < p)) print $1 }' $file >> te_fimo_sequences.txt
done
sort -u te_fimo_sequences.txt > tmp; mv tmp te_fimo_sequences.txt
seqkit grep -f te_fimo_sequences.txt  ../../data/genome/ucsc_hg38_reps.fa > te_fimo_sequences.fa

cat ../../data/cis-reg/kzfp_candidate_promoters.fa > kzfp_fimo_sequences.fa

# Run FIMO to identify expressed iPSC TF motifs against TE and KZFP sequences
fimo \
    -oc ../../data/cis-reg/te-fimo \
    -bfile genome_background.bg \
    ../../data/cis-reg/ipsc_motifs.meme \
    te_fimo_sequences.fa
fimo \
    -oc ../../data/cis-reg/kzfp-fimo \
    -bfile genome_background.bg \
    ../../data/cis-reg/ipsc_motifs.meme \
    kzfp_fimo_sequences.fa

# rm genome_background.bg
# rm te_fimo_sequences.{fa,txt}
# rm kzfp_fimo_sequences.{fa,txt}
# rm motif_dict.txt search_names.txt mincount_names.txt omit_kzfps.txt meme_names.txt
