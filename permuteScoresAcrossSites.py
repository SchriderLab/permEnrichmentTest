import sys, os, random, math

class Chrom:
    def __init__(self, name, sites):
        self.name = name
        self.sites = sites


def readClrScoreFile(clrScoreFileDir, chromNameSuffix):
    chroms = []
    clrScores = []
    clrScoreFiles = os.listdir(clrScoreFileDir)
    for clrScoreFile in clrScoreFiles:

        chromName = clrScoreFile.split(chromNameSuffix)[0]
        with open(clrScoreFileDir + "/" + clrScoreFile, "rt") as csf:

            first = True
            sites = []
            for line in csf:
                if first:
                    header = line.strip()
                    first = False
                else:
                    line = line.strip().split()
                    sites.append(line[0] + "\t" + line[1])
                    score = line[2] + "\t" + line[3]
                    clrScores.append(score)

        chroms.append(Chrom(chromName, sites))

    return header, chroms, clrScores


def permScores(clrScores, chunkSize=None):
    permedClrScores = []

    startIndex = random.randint(0, len(clrScores)-1)

    # segment the genome into chunks and randomly shuffle them
    if chunkSize:
        bins = []
        numBins = int(math.ceil(len(clrScores) / chunkSize))
        for binIndex in range(numBins):
            bins.append([])

        for i in range(len(clrScores)):
            binIndex = int(i / chunkSize)
            bins[binIndex].append(clrScores[i])
        permedBins = random.sample(bins, k=len(bins))

        for binIndex in range(numBins):
            for score in permedBins[binIndex]:
                permedClrScores.append(score)
        assert len(clrScores) == len(permedClrScores)
    # giant concatenated chromosome with a random starting point
    else:
        for i in range(len(clrScores)):
            currIndex = (i+startIndex) % len(clrScores)
            permedClrScores.append(clrScores[currIndex])

    return permedClrScores


def main():

    clrScoreFileDir, chromNameSuffix, numPerms, chunkSize, outDir = sys.argv[1:]
    numPerms = int(numPerms)
    chunkSize = int(chunkSize)

    if chunkSize < 0:
        sys.exit("chunkSize must be a non-negative intger. (See README)")

    header, chroms, clrScores = readClrScoreFile(clrScoreFileDir, chromNameSuffix)

    sys.stderr.write("starting permutations\n")
    for perm in range(numPerms):
        permedClrScores = permScores(clrScores, chunkSize=int(chunkSize))

        i = 0
        for chrom in random.sample(chroms, k=len(chroms)):
            permOutFile = f"{outDir}/{chrom.name}_perm_{perm}.txt"
            with open(permOutFile, 'wt') as pof:
                pof.write(header + "\n")
                for site in chrom.sites:
                    pof.write(f"{site}\t{permedClrScores[i]}\n")
                    i += 1
        sys.stderr.write(f"done {perm} of {numPerms} permutations--------------\r")
    sys.stderr.write("\ndone\n")

if __name__ == "__main__":
    main()
