# Open condensed GWAS file
gwasFile = 'Netblotch.cgwas'
try:
    gwasFh = open(gwasFile, 'r')
except IOError:
    print('There was an error opening the condensed GWAS input file')
    exit(1)

gffFile = 'Hv_IBSC_PGSB_v2p37_genes.txt'
try:
    gffFh = open(gffFile, 'r')
except IOError:
    print('There was an error opening the gff (.txt) file')
    exit(1)

# Save the ranges of interest to dictionary
gwasRanges = {} # {chrom: [begin, end]}
keys = []
for line in gwasFh:
    chrom = line.split()[1]
    if chrom not in keys:
        keys.append(chrom)
        gwasRanges[chrom] = []
    lineSplit = line.split()
    del lineSplit[4], lineSplit[1]
    gwasRanges[chrom].append(lineSplit)
gwasFh.close()

outFile = gwasFile.replace('.cgwas', '.gff3.out')
try:
    outFh = open(outFile, 'w')
except IOError:
    print('There was an error opening the output file')
    exit(1)

# Find and save all gff ids for ranges that lie within target ranges
snpGeneOverlaps = {}
snpKeys = []
for line in gffFh:
    current = line.split()
    if current[0] in keys:
        [chrom, gffBegin, gffEnd, ID] = current
        [gffBegin, gffEnd] = [int(gffBegin), int(gffEnd)]
        # beginning or end lies between the current gwas range, save it to file
        for entry in gwasRanges[chrom]:
            [snp, rangeBegin, rangeEnd] = [entry[0], int(entry[1]), int(entry[2])]
            if (gffBegin >= rangeBegin and gffBegin <= rangeEnd) or (gffEnd >= rangeBegin and gffEnd <= rangeEnd):
                if snp not in snpKeys:
                    snpKeys.append(snp)
                    snpGeneOverlaps[snp] = []
                snpGeneOverlaps{}
                print(entry, current)
                # out.write(line)

gffFh.close()
outFh.close()
