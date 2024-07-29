#!/usr/bin/env python3

import re
import pybedtools as pb
from scipy import stats, signal
from collections import defaultdict

"""
This script extracts ZNF clusters from bed/gff files of gene coordinates and prints to bedfile
format.
"""



class ZFGenome:

    def __init__(self, zf_bedfile, genomefile):
        self.genomefile = genomefile
        self._idx = 0
        self._chroms = set()
        self._clust_idx = 0
        zf_bed = pb.BedTool(zf_bedfile)
        self.zf_bed = zf_bed.each(self._rename)

    def _rename(self, feature):
        feature.name = f'{feature.name}_{self._idx}'
        self._idx += 1
        return feature

    def _parse_cluster(self, feature):
        if feature.chrom not in self._chroms:
            self._chroms.add(feature.chrom)
            self._clust_idx = 0
        feature.name = f'{feature.chrom}_{self._clust_idx}'
        self._clust_idx += 1
        feature.strand = '.'
        return feature

    def extract_zf_clusters(self, maxdist=5e4, margin=1e4):
        """Use intra-gene distance to assign KZFPs to clusters.
        
        Arguments:
            maxdist: maximum distance between genes to merge into cluster.
            margin: size of margins around each cluster
        """

        maxdist = int(maxdist)
        margin = int(margin)
        # HACK: redundant 'counts' needed to insure BED intervals have
        merged_bed = self.zf_bed.merge(d=maxdist, 
                                       c=(4, 4, 4, 4), 
                                       o=('count', 'count', 'count', 'collapse')) \
                                .slop(b=margin, g=self.genomefile)
        self.merged_bed = merged_bed.each(self._parse_cluster)
            
    def __repr__(self):
        return str(self.merged_bed)

if __name__ == '__main__':
    with open('../../data/parsed_metazoans.out') as infile:
        for line in infile:
            species = line.split('\t')[0]
            print(species)
            assembly_type = line.strip().split('\t')[-1]
            if assembly_type != 'Chromosome':
                continue
            try:
                zf_genome = ZFGenome(f'/Users/jonwells/Projects/feschottelab/metazoan-znfs/data/beds/{species}_znfs.bed',
                                     f'/Users/jonwells/Projects/feschottelab/metazoan-znfs/data/genomes/{species}.genome')
                zf_genome.extract_zf_clusters(maxdist=1e5)
                with open(f'../../data/beds/{species}_zf_clusters.bed', 'w') as outfile:
                    outfile.write(str(zf_genome))
            except:
                print(f'WARNING: issue parsing {species}')



