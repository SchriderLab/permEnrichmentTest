#!/usr/bin/env python3

import os, sys
import argparse
import numpy as np

def readGeneCoords(geneCoordFileName):
    genes = []
    chromNames = {'AaegL5_1': 'AaegL5_1', 'AaegL5_2': 'AaegL5_2', 'AaegL5_3': 'AaegL5_3', 'NC_035107.1': 'AaegL5_1', 'NC_035108.1': 'AaegL5_2', 'NC_035109.1': 'AaegL5_3'}
    with open(geneCoordFileName, 'rt') as f:
        for line in f:
            if not (line.startswith("#") or line.startswith("NIGP")):
                line = line.strip().split("\t")

                if line[0].startswith("Aaeg") or line[0].startswith("NC_"):
                    chrom, start, end, info = line
                    try:
                        geneName = info.split(";Name=")[1].split(";")[0]
                    except Exception:
                        geneName = info.split("ID=")[1].split(";")[0]
                else:
                    geneId, geneName, pub, chrom, start, end = line

                if chrom != ".":
                    chrom = chromNames[chrom]
                    genes.append((chrom, int(start), int(end), geneName))
    return genes

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

    for chrom, start, end, name in genes:
        for pos in range(start-bufferDist, end+bufferDist+1):
            if (chrom, pos) in siteToGenes:
                siteToGenes[(chrom, pos)][name] = 1

    return siteToGenes


#TODO: need output args, and also our input GO def file arg
parser = argparse.ArgumentParser()
parser.add_argument("--real_scores", help="path to file with real input scores", required=True)
parser.add_argument("--perm_scores", help="path to directory with permuted input scores", required=True)
parser.add_argument("--perm_delimiter", help="path to directory with permuted input scores", default="_perm_")
parser.add_argument("--outlier_percentile", 
                    help="int or float that specifies what score percentile to use as the cutoff for outliers (e.g. 97.5)",
                    type=float,
                    required=True)
parser.add_argument("--n_perms", help="number of permutations", type=int, required=True)
parser.add_argument("--annot_file", help="path to file with gene names and coordinates", required=True)
parser.add_argument("--gene_buffer_dist", help="sites within this distance of a gene are counted as within that gene", type=int, default=0)

args = parser.parse_args()


sys.stderr.write("reading in scores\n")
scoresAndSites = readScoresAndSites(args.real_scores)

sites = [x[0] for x in scoresAndSites]
chroms = [x[0] for x in sites]
chroms = set(chroms)

sys.stderr.write("reading gene coords\n")
# read in gene info from our gene list
genes = readGeneCoords(args.annot_file)


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
#outlierSites = [x[0] for x in scoresAndSites]
outlierSites = [x[0] for x in scoresAndSites if x[1] > cutoff]

#realOutlierGenes = {}
#for site in outlierSites:
#    for gene in siteToGenes[site]:
#        realOutlierGenes[gene] = 1
#print("outlier genes in real data: ")
#for gene in realOutlierGenes:
#    print(gene)

def posToWin(pos):
    return pos - (pos % 1000000)

realOutlierGeneWins = {}
for site in outlierSites:
    if len(siteToGenes[site]) > 0:
        chrom, pos = site
        win = (chrom, posToWin(pos))
        realOutlierGeneWins[win] = 1
print("outlier gene windows in real data: ")
for win in realOutlierGeneWins:
    print(win)

#realOutlierCount = len(realOutlierGenes)
realOutlierCount = len(realOutlierGeneWins)

print("total number of outlier gene wins in real data: ", realOutlierCount)

# count the number of distinct genes each term occurs in in our outlier set
# (note that this includes zero counts)
# first initialize
sys.stderr.write("getting gene wins with at least one outlier score in real data\n")


sys.stderr.write("examining permuted scores\n")
# get outlier genes based on our permuted scores
allPermOutlierCounts = []
for i in range(args.n_perms):
#for i in range(100):

    outlierSites = []
    for chrom in chroms:
        permFileNameForChrom = f"{chrom}{args.perm_delimiter}{i}.txt"
        outlierSites += [x[0] for x in readScoresAndSites(args.perm_scores+"/"+permFileNameForChrom) if x[1] > cutoff]

    #permOutlierGenes = {}
    #for site in outlierSites:
    #    for gene in siteToGenes[site]:
    #        permOutlierGenes[gene] = 1
    #allPermOutlierCounts.append(len(permOutlierGenes))

    permOutlierGeneWins = {}
    for site in outlierSites:
        if len(siteToGenes[site]) > 0:
            chrom, pos = site
            win = (chrom, posToWin(pos))
            permOutlierGeneWins[win] = 1
    allPermOutlierCounts.append(len(permOutlierGeneWins))

    #sys.stderr.write(f"done getting outlier genes and terms for {i} of {args.n_perms} permuted score sets----------------\r")
    sys.stderr.write(f"done getting outlier gene wins and terms for {i} of {args.n_perms} permuted score sets----------------\r")
sys.stderr.write("\nall done!\n")

# get our one-sided p-val
numPermsGeReal = len([x for x in allPermOutlierCounts if x >= realOutlierCount])
pval = numPermsGeReal / len(allPermOutlierCounts)

# get our fold-enrichment
meanPerm = np.mean(allPermOutlierCounts)
if meanPerm == 0:
    enrichment = 0
else:
    enrichment = realOutlierCount / meanPerm

print(f"\nreal intersect: {realOutlierCount}; mean permuted intersect: {meanPerm}; enrichment: {enrichment}; p-value: {pval}")
