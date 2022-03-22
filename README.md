# permEnrichmentTest
This respository is where we will host some scripts for doing permutation-based enrichment tests.

## Permuting scores (e.g. from a sweep scan) across a genome while accounting for linkage and variable gene lengths
If we wish to do a permutation test for enrichment that accounts for the autocorrelation of scores, for vairation in gene lengths, and also the physical clustering of genes along the chromosome, we can do this by constructing a concatenated cicrularized chromosome (following https://arxiv.org/abs/2004.03224), and then simply shifting our scores along this circularized chromosome by some random value (permutation strategy 1).

If we have a small number of chromosomes and/or relatively wide signal peaks, there may be some concern that this may not give us enough randomness (i.e. roughly the same arrangement of scores in some of the permuted datasets as the real data). Therefore, if desired we can go a step further by partioning the chromosome into large chunks, and then shuffling the order of these chunks to form our permuted datasets (permutation strategy 2). Both of these permutation strategies are implemented in `permuteScoresAcrossSites.py`, which can be run as follows:

```python permuteScoresAcrossSites.py scores/ chromNameDelimiter numPerms chunkSize permutedScores/```

where `scores/` is the path to a directory that contains our real results (format discussed below), `chromNameDelimiter` is explained just below, `numPerms` is the number of permutations we wish to perform, `chunkSize` is the number of consecutive scores assigned to each chunk when using permutation strategy 2 described above, and `permutedScores` is the directory where we will write all of our output. If `chunkSize` is set to 0, we permute using strategy 1, otherwise we use strategy 2.

This script assumes that our input directory contains one score file per chormosome, and that the name of each score file begins with the name of the chromosome. `chromNameDelimiter` is a string specifying where in this file name the chromosome name ends. For example, if our input file is `AaegL5_1_merged.031622.csv`, the chromosome name is `AaegL5_1`, so we will supply `_merged` as the second argument to this script (i.e. as `chromNameDelimiter`). The output directory will then contain `numPerms` files for each chromosome, consisting of permuted scores assigned to the same positions as specified in our original `scores/` directory.

Currently, `scores` takes an input format that is identical to the output from SweepFinder, but with an additional first column that specifies the chromosome name. Thus, this is a tab-delimited file with the following four fields: chromosome name ("chrom" in the header line), the tested site ("location"), SweepFinder's likelihood ratio ("LR"), and SweepFinder's alpha estimate ("alpha").
