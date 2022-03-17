#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:35:18 2022

@author: beccalove
"""

import argparse
import numpy as np
import pandas as pd
import random

import allel

from permuteScoresAcrossSites import permScores

def findGenesWithOutliers(scores, genes, cutoff):
    
    outlier_sites = scores.loc[scores["LR"] > cutoff]
    
    chroms = scores["chrom"].unique()
    
    genes_with_outlier = []

    for chrom in chroms:

        outlier_chunk = outlier_sites.loc[outlier_sites["chrom"] == chrom]

        genes_chunk = genes.loc[genes["seqid"] == chrom]

        _, overlapped_genes_bool =\
        allel.SortedIndex(outlier_chunk["location"]).locate_intersection_ranges(
            genes_chunk["start"], genes_chunk["end"])

        genes_with_outlier.extend(genes_chunk.loc[overlapped_genes_bool, "name"].values)

    return genes_with_outlier

parser = argparse.ArgumentParser()
parser.add_argument("--scores", help="path to input scores")
parser.add_argument("--gff3", help="path to genes as defined in .gff3 file")
parser.add_argument("--outlier_percentile", 
                    help="int or float that specifies what score percentile to use as the cutoff for outliers (ex. 95)",
                    type=float)
parser.add_argument("--n_perms", help="number of permutations", type=int)

args = parser.parse_args()

##def main():

    ##read in and munge data
scores = pd.read_table(args.scores)

gff3 = pd.read_table(args.gff3, sep="\t", comment="#", header=None)

gff3.columns = ["seqid", "source", "type", "start", "end", "score", "strand",
               "phase", "attributes"]

gff3["name"] = gff3["attributes"].str.split(";", 
    expand=True)[0].str.lstrip("ID=")

genes = gff3.loc[gff3["type"] == "protein_coding_gene"]

cutoff = np.percentile(scores["LR"].values, args.outlier_percentile)

##find genes with outliers in original data, then permute
originals = findGenesWithOutliers(scores, genes, cutoff)

permuted = []

for i in range(args.n_perms):
    
    new_scores = scores.copy()
    
    new_scores["LR"] = permScores(scores["LR"].values)
    
    permuted.append(findGenesWithOutliers(new_scores, genes, cutoff))
        
##if __name__ == "__main__":
##    main()