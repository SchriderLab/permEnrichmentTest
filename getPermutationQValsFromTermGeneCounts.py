import sys,os

realCountFile, permCountDir = sys.argv[1:]

def readTermSigCounts(countFileName, counts, recordDefs=False):
    if recordDefs:
        termDefs = {}

    with open(countFileName) as intersectFile:
        for line in intersectFile:
            term, termDef, numTermSigGenes, termSigGenes = line.rstrip("\n").split("\t")
            count = int(numTermSigGenes)
            if not term in counts:
                counts[term] = 0
            counts[term] += count

            if recordDefs:
                termDefs[term] = termDef

    if recordDefs:
        return termDefs

permCounts = []

realCounts = {}
sys.stderr.write("reading real term counts\n")
termNames = readTermSigCounts(realCountFile, realCounts, recordDefs=True)

sys.stderr.write("starting to read perm term counts\n")
permFileNames = os.listdir(permCountDir)
for i in range(len(permFileNames)):
    currPermCounts = {}

    for term in realCounts:
        currPermCounts[term] = 0

    readTermSigCounts(permCountDir + "/" + permFileNames[i], currPermCounts)
    permCounts.append(currPermCounts)
    sys.stderr.write(f"done {i} of {len(permFileNames)} perms--------\r")
sys.stderr.write("\ndone\n")

outLs = []
totalCounts = {}
outLineH = {}
sys.stderr.write("calculating p-values\n")
i = 0
for term in realCounts:
    pCount = 0
    totalCount = 0
    permutedIntersectSum = 0
    for permCountH in permCounts:
        if term in permCounts:
            permCount = permCounts[term]
        else:
            permCount = 0
        permutedIntersectSum += permCount
        if permCount >= realCounts[term]:
            pCount += 1
        totalCount += 1
    totalCounts[totalCount] = 1

    meanIntersect = permutedIntersectSum/float(totalCount)
    if meanIntersect > 0:
        enrichment = realCounts[term] / meanIntersect
    else:
        enrichment = float("inf")
    if pCount == 0:
        pValStr = "<%s" %(1.0/totalCount)
    else:
        pValStr = str(pCount / float(totalCount))
    if not pCount in outLineH:
        outLineH[pCount] = []

    outLineH[pCount].append((enrichment, "%s: %s; real intersect: %s; mean permuted intersect: %s; enrichment: %s; p-value: %s" %(term, termNames.get(term, term), realCounts[term], meanIntersect, enrichment, pValStr)))

    sys.stderr.write(f"done {i} of {len(realCounts)} terms-----------------\r")
    i += 1

assert len(totalCounts) == 1
totalCount = list(totalCounts.keys())[0]
sys.stderr.write(f"\ndone\n")


positiveCount = len(realCounts)
minQVal=1.0
minNonZeroQVal=1.0
sys.stderr.write("calculating q-values\n")
for pCount in sorted(outLineH, reverse=True):
    pVal = pCount/float(totalCount)
    fdr = (pVal*len(realCounts))/positiveCount
    if fdr < minQVal:
        if fdr != 0:
            minNonZeroQVal = fdr
        minQVal = fdr
    if minQVal > 0:
        qValStr = str(minQVal)
    else:
        qVal = ((1/float(totalCount))*len(realCounts))/positiveCount
        if qVal < minNonZeroQVal:
            qValStr = "<%s" %(qVal)
        else:
            qValStr = "<%s" %(minNonZeroQVal)
    outLines = []
    for enrichment, outLine in sorted(outLineH[pCount]):
        print(outLine + "; q-value: %s" %(qValStr))
    positiveCount -= len(outLineH[pCount])

sys.stderr.write("all done!\n")
