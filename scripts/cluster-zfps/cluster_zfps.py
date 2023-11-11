#!/usr/bin/env python3

import re
import pybedtools as pb
from scipy import stats, signal
from collections import defaultdict

def load_chr_lengths(genomefile):
    with open(genomefile) as infile:
        chr_lengths = {l.split()[0]: (0, int(l.split()[1])) for l in infile}
    return chr_lengths

def load_chom_accessions(accessions_file):
    with open(accessions_file) as infile:
        acc_to_chr = {l.split()[0]: l.split()[1] for l in infile}
    return acc_to_chr    

def load_zfp_coordinates(zfp_transcript_gff):
    zfp_gff = pb.BedTool(zfp_transcript_gff)
    zfp_coords = defaultdict(dict)
    for feature in zfp_gff:
        midpoint = int(feature.start + (feature.end - feature.start)/2)
        zfp_coords[feature.chrom][feature.name.strip('rna-')] = (feature.start, midpoint, feature.end)
    return zfp_coords

def extract_zfp_clusters(zfp_coords, chr_lengths, min_cluster_size=3, acc_to_chr=None, bandwidth=0.025):
    """Use kernel density estimation to assign KZFPs to clusters.
    
    Arguments:
        zfp_coords - dictionary containing start, midpoint and end of each zfp_transcript
        chr_lengths - dictionary containing lengths for each chromosome
        bandwidth - KDE bandwith, i.e. cluster granularity
    Returns:
        zfp_clusters - dictionary mapping from cluster to list of constituent genes
    """
    zfp_clusters = {}
    for chr_accession in zfp_coords.keys():
        if len(zfp_coords[chr_accession]) < min_cluster_size:
            continue
        zfp_kde = stats.gaussian_kde([v[1] for v in zfp_coords[chr_accession].values()], 
                                      bw_method=bandwidth)
        start, stop, step = *chr_lengths[chr_accession], 10000
        zfp_pdf = zfp_kde.pdf([*range(start, stop, step)])
        zfp_pdf_minima = [start] + list(step*signal.find_peaks(-zfp_pdf)[0]) + [stop]
        
        # Assign genes to clusters using gend KDE minima as boundaries
        wide_cluster_boundaries = [(zfp_pdf_minima[i], zfp_pdf_minima[i+1]) for i in range(len(zfp_pdf_minima)-1)]
        chrom_zfp_clusters = []
        curr_cluster = []
        i = 0
        for gene, mp in zfp_coords[chr_accession].items():
            mp = mp[1]
            start, stop = wide_cluster_boundaries[i]
            if mp > stop:
                # Increment wide cluster and reset current cluster list to empty
                i += 1
                chrom_zfp_clusters.append(curr_cluster)
                curr_cluster = []
            curr_cluster.append(gene)
        chrom_zfp_clusters.append(curr_cluster)
        
        # Filter out clusters not meeting minimum size threshold and label: chrom.cluster 
        for i, cluster in enumerate([c for c in chrom_zfp_clusters if len(c) >= min_cluster_size], 1):
            if acc_to_chr:
                zfp_clusters[f'{acc_to_chr[chr_accession]}.{i}'] = cluster
            else:
                zfp_clusters[f'{chr_accession}.{i}'] = cluster
    
    return zfp_clusters

def extract_cluster_boundaries(clusters, zfp_coords, chr_lengths, margin=2e5):
    """Extact genomic coordinates of each cluster plus flanking (margin) sequence. 
    
    Arguments:
        clusters - dictionary containing genes in each cluster
        zfp_coords - dictionary containing start, midpoint and end of each zfp_transcript
        margin - slop from the edges of first and last genes in cluster
    Returns:
        cluster_boundaries - dictionary mapping start and end coordinates of clusters
    """
    cluster_boundaries = {}
    flat_zfp_coords = {}
    for chrom, zfpdict in zfp_coords.items():
        for zfp, (start, midpoint, end) in zfpdict.items():
            flat_zfp_coords[zfp] = (chrom, start, end)
    for key, cluster in clusters.items():
        chrom = flat_zfp_coords[cluster[0]][0]
        coords = [flat_zfp_coords[zfp][1:] for zfp in cluster]
        coords = [i for j in coords for i in j]
        # Set cluster boundaries as egde of first/last gene plus the margin, unless this exceeds chromosome boundaries.
        cluster_boundaries[key] = (chrom, 
                                   int(max((min(coords)-margin, 0))),
                                   int(min((max(coords)+margin, chr_lengths[chrom][1]))))
    return cluster_boundaries

def cluster_to_bed(cluster_boundaries, bedfile):
    with open(bedfile, 'w') as outfile:
        for clustername, (chrom, start, end) in cluster_boundaries.items():
            outfile.write(f'{chrom}\t{start}\t{end}\t{clustername}\n')
        

if __name__ == '__main__':
    chr_lengths = load_chr_lengths('/Users/jonwells/Genomes/Mollusca/GCF_006345805.1_ASM634580v1.genome')
    # acc_to_chr = load_chom_accessions('../../data/GCF_009914755.1.chrom_accessions.txt')
    # chr_to_acc = {v: k for (k, v) in acc_to_chr.items()}
    zfp_coords = load_zfp_coordinates('../search-kzfps/octopus_transcripts_only.gff')
    zfp_clusters = extract_zfp_clusters(zfp_coords, chr_lengths)
    cb = extract_cluster_boundaries(zfp_clusters, zfp_coords, chr_lengths)
    cluster_to_bed(cb, './octopus_test.bed')
    # for zfp in zfp_clusters['19.10']:
    #     start, mp, stop = zfp_coords[chr_to_acc['19']][zfp]
    #     print(f'{zfp}\t{start}\t{mp}\t{stop}')


