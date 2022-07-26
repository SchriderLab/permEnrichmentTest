#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:35:18 2022

@author: beccalove
"""

import os, sys
import argparse
import numpy as np
import pandas as pd
import random

#import allel

#def findGenesWithOutliers(scores, genes, cutoff):
#    
#    outlier_sites = scores.loc[scores["LR"] > cutoff]
#    
#    chroms = scores["chrom"].unique()
#    
#    genes_with_outlier = []
#
#    for chrom in chroms:
#
#        outlier_chunk = outlier_sites.loc[outlier_sites["chrom"] == chrom]
#
#        genes_chunk = genes.loc[genes["seqid"] == chrom]
#
#        _, overlapped_genes_bool =\
#        allel.SortedIndex(outlier_chunk["location"]).locate_intersection_ranges(
#            genes_chunk["start"], genes_chunk["end"])
#
#        genes_with_outlier.extend(genes_chunk.loc[overlapped_genes_bool, "name"].values)
#
#    return genes_with_outlier


def readScoresAndSites(scoreFileName):
    scoresAndSites = []
    first = True

    with open(scoreFileName, 'rt') as sf:
        for line in sf:
            if first:
                first = False
            else:
                chrom, pos, clr, alpha = line.strip().split("\t")[:4]
                pos = int(pos)
                clr = float(clr)
                site  = (chrom, pos)
                scoresAndSites.append((site, clr))

    return scoresAndSites


def getSiteToGenes(sites, genes, bufferDist):
    siteToGenes = {}
    for chrom, pos in sites:
        siteToGenes[(chrom, pos)] = {}

    for index, gene in genes.iterrows():
        chrom = gene['seqid']
        for pos in range(gene['start']-bufferDist, gene['end']+bufferDist+1):
            if (chrom, pos) in siteToGenes:
                siteToGenes[(chrom, pos)][gene['name']] = 1

    return siteToGenes


def readAnnotFile(annotFileName):
    first = True
    termDefs = {}
    geneToTerms = {}

    with open(annotFileName, 'rt') as af:
        for line in af:
            if first:
                first = False
            else:
                geneName, termId, termDesc = line.strip().split("\t")
                termDefs[termId] = termDesc
                if not geneName in geneToTerms:
                    geneToTerms[geneName] = []
                geneToTerms[geneName].append(termId)

    return geneToTerms, termDefs

#TODO: need output args, and also our input GO def file arg
parser = argparse.ArgumentParser()
parser.add_argument("--real_scores", help="path to file with real input scores", required=True)
parser.add_argument("--perm_scores", help="path to directory with permuted input scores", required=True)
parser.add_argument("--perm_delimiter", help="path to directory with permuted input scores", default="_perm_")
parser.add_argument("--gff3", help="path to genes as defined in .gff3 file", required=True)
parser.add_argument("--outlier_percentile", 
                    help="int or float that specifies what score percentile to use as the cutoff for outliers (e.g. 97.5)",
                    type=float,
                    required=True)
parser.add_argument("--n_perms", help="number of permutations", type=int, required=True)
parser.add_argument("--annot_file", help="path to tab-delimited file with gene names, annotated terms, and term descriptions", required=True)
parser.add_argument("--real_out_file", help="path where annotation counts for the real data will be written", required=True)
parser.add_argument("--perm_out_dir", help="directory where annotation counts for permuted data will be written", required=True)
parser.add_argument("--gene_buffer_dist", help="sites within this distance of a gene are counted as within that gene", type=int, default=0)

args = parser.parse_args()


sys.stderr.write("reading in scores\n")
scoresAndSites = readScoresAndSites(args.real_scores)

sites = [x[0] for x in scoresAndSites]
chroms = [x[0] for x in sites]
chroms = set(chroms)

sys.stderr.write("reading gene coords from gff\n")
# read in gene info from our gff3 file
gff3 = pd.read_csv(args.gff3, sep="\t", comment="#", header=None)
gff3.columns = ["seqid", "source", "type", "start", "end", "score", "strand",
               "phase", "attributes"]
gff3["name"] = gff3["attributes"].str.split(";", 
    expand=True)[0].str.lstrip("ID=")
genes = gff3.loc[gff3["type"] == "protein_coding_gene"]


sys.stderr.write(f"mapping tested sites to genes with a buffer of {args.gene_buffer_dist} bp on either side\n")
# map sites to genes, and get all tested genes
siteToGenes = getSiteToGenes(sites, genes, args.gene_buffer_dist)

allTestedGenes = {}
for key in siteToGenes.keys():
    for gene in siteToGenes[key]:
        allTestedGenes[gene] = 1

print("total tested genes: ", len(allTestedGenes))

cutoff = np.percentile([x[1] for x in scoresAndSites], args.outlier_percentile)
print("score cutoff: ", cutoff)

sys.stderr.write("getting outlier scores from real data\n")
# get outlier genes based on our real scores
outlierSites = [x[0] for x in scoresAndSites if x[1] > cutoff]
outlierGenes = {}
for site in outlierSites:
    for gene in siteToGenes[site]:
        outlierGenes[gene] = 1

print("outlier genes in real data: ", len(outlierGenes))

sys.stderr.write("reading term annotations\n")
# read annotations
geneToTerms, termDefs = readAnnotFile(args.annot_file)

# count the number of distinct genes each term occurs in in our outlier set
# (note that this includes zero counts)
# first initialize
sys.stderr.write("getting genes with at least one outlier score in real data\n")
realTermSigGenes = {}
for gene in allTestedGenes:
    if gene in geneToTerms:
        for term in geneToTerms[gene]:
            realTermSigGenes[term] = []

# now record the outlier genes
for gene in outlierGenes:
    if gene in geneToTerms:
        for term in geneToTerms[gene]:
            realTermSigGenes[term].append(gene)


# write them out!
with open(args.real_out_file, 'wt') as of:
    for term in realTermSigGenes:
        of.write(f"{term}\t{termDefs[term]}\t{len(realTermSigGenes[term])}\t{','.join(realTermSigGenes[term])}\n")

sys.stderr.write("examining permuted scores\n")
# get outlier genes based on our permuted scores
for i in range(args.n_perms):

    outlierSites = []
    for chrom in chroms:
        permFileNameForChrom = f"{chrom}{args.perm_delimiter}{i}.txt"
        outlierSites += [x[0] for x in readScoresAndSites(args.perm_scores+"/"+permFileNameForChrom) if x[1] > cutoff]

    outlierGenes = {}
    for site in outlierSites:
        for gene in siteToGenes[site]:
            outlierGenes[gene] = 1

    # first initialize
    permTermSigGenes = {}
    for gene in allTestedGenes:
        if gene in geneToTerms:
            for term in geneToTerms[gene]:
                permTermSigGenes[term] = []

    # now record our outlier genes
    for gene in outlierGenes:
        if gene in geneToTerms:
            for term in geneToTerms[gene]:
                permTermSigGenes[term].append(gene)

    # write them out!
    outFilePath = f"{args.perm_out_dir}/terms{args.perm_delimiter}{i}.txt"
    with open(outFilePath, 'wt') as of:
        for term in permTermSigGenes:
            of.write(f"{term}\t{termDefs[term]}\t{len(permTermSigGenes[term])}\t{','.join(permTermSigGenes[term])}\n")
    sys.stderr.write(f"done getting outlier genes and terms for {i} of {args.n_perms} permuted score sets----------------\r")
sys.stderr.write("\nall done!\n")
