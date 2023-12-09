These scripts were used to generate an approximation to the human KZFP-TE gene regulatory network,
as found in hiPSCs. This was constructed according to the following protocol:

1. map_ipsc_reads.sh 
   Reads from a human induced pluripotent stem cell line (HipSci database: HPSI0114i-bezi_1;
   ERR947017) were aligned to the human genome (T2T-CHM13v2.0; GCF_009914755.1) using STAR v2.7.11a.
   Reads were counted across both genes and TEs using TEcounts.

2. extract_te_cres.sh
   To identify candidate promoters an enhancers for KZFPs, and to identify TE-associated TRIM28/KAP1
   binding sites within the vicinity of KZFPs. This enables edges to be drawn from KZFPs to other
   KZFPs.

3. extract_putative_tfs.sh
   A list of TFs active in hiPSCs was selected based on a minimum read count threshold (see script),
   and the DNA-binding motifs associated with these TFs were collected from the JASPAR 2022
   non-redundant core database. These motifs were scanned against KZFP candidate promoter sequences
   and any TEs significantly targeted by at least one KZFP (see script for details).

4. extract_zfp_targets.sh
   In step two, TRIM28 peaks lying within the vicinity of KZFPs are associated with specific TEs.
   However, many peaks are KZFP-associated but do not overlap a particular TE insertion. We
   therefore assess KZFP targeting to these peaks by scanning for binding motifs. This enables
   identification of direct KZFP -| KZFP interactions.

5. parse_cisbp.py
   Helper script that is used within step four to convert from CisBP PWM format to MEME motifs.

6. cres_to_network.py
   This script provides functions that can be used to parse and combine the output of those 
   previously described, and ultimately generates a simple adjacency list that can be used to model
   a directed GRN.
