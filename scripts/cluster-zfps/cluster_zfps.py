#!/usr/bin/env python3

import re
import pybedtools as pb
from scipy import stats
from collections import defaultdict

"""
This script extracts ZNF clusters from bed/gff files of gene coordinates and prints to bedfile
format.
"""


class ZFGenome:

    def __init__(self, zf_bedfile, genomefile, use_midpoint=True):
        self.genomefile = genomefile
        self._idx = 0
        self._chroms = set()
        self._clust_idx = 0
        self.zf_bed = pb.BedTool(zf_bedfile)
        self.zf_bed = self.zf_bed.each(self._rename)
        if use_midpoint:
            self.zf_bed = self.zf_bed.each(self._extract_midpoint)

    def _rename(self, feature):
        feature.name = f'{feature.name}_{self._idx}'
        self._idx += 1
        return feature

    def _extract_midpoint(self, feature):
        assert feature.start <= feature.stop
        midpoint = feature.start + (feature.stop-feature.start)//2
        feature.start = midpoint
        feature.stop = midpoint + 1
        return feature

    def _parse_cluster(self, feature):
        if feature.chrom not in self._chroms:
            self._chroms.add(feature.chrom)
            self._clust_idx = 0
        feature.name = f'{feature.chrom}_{self._clust_idx}'
        self._clust_idx += 1
        feature.strand = '.'
        return feature

    def extract_zf_clusters(self, maxdist, margin=1e4):
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
                                       o=('count', 'count', 'count', 'collapse'))
        merged_bed = merged_bed.each(self._parse_cluster)
        self.merged_bed = merged_bed.slop(b=margin, g=self.genomefile)
            
    def __repr__(self):
        return str(self.merged_bed)

if __name__ == '__main__':

    max_inter_znf_dist = 118261 # 75th percentile
    use_midpoint = False

    with open('../../data/parsed_metazoans.out') as infile:
        for line in infile:
            species = line.split('\t')[0]
            assembly_type = line.strip().split('\t')[-1]
            if assembly_type not in ['Chromosome']:
                continue
            try:
                print(species)
                zf_genome = ZFGenome(f'/Users/jonwells/Projects/feschottelab/metazoan-znfs/data/beds/{species}_znfs.bed',
                                     f'/Users/jonwells/Projects/feschottelab/metazoan-znfs/data/genomes/{species}.genome',
                                     use_midpoint=use_midpoint)
                zf_genome.extract_zf_clusters(maxdist=max_inter_znf_dist)
                with open(f'../../data/beds/{species}_zf_clusters.bed', 'w') as outfile:
                    outfile.write(str(zf_genome))
            except:
                print(f'WARNING: issue parsing {species}')



