# Open condensed GWAS file
gwasFile = 'water_absorption.cgwas'
try:
    gwasFh = open(gwasFile, 'r')
except IOError:
    print('There was an error opening the condensed GWAS input file')
    exit(1)

gffFile = 'Triticum_aestivum.TGACv1.38.txt'
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
    gwasRanges[chrom].append(line.split()[2:4])
gwasFh.close()

outFile = gwasFile.replace('.cgwas', '.gff3.out')
try:
    outFh = open(outFile, 'w')
except IOError:
    print('There was an error opening the output file')
    exit(1)

# Find and save all gff ids for ranges that lie within target ranges
for line in gffFh:
    current = line.split()
    if current[0] in keys:
        [chrom, gffBegin, gffEnd, ID] = current
        [gffBegin, gffEnd] = [int(gffBegin), int(gffEnd)]
        # beginning or end lies between the current gwas range, save it to file
        for entry in gwasRanges[chrom]:
            [rangeBegin, rangeEnd] = [int(entry[0]), int(entry[1])]
            if (gffBegin >= rangeBegin and gffBegin <= rangeEnd) or (gffEnd >= rangeBegin and gffEnd <= rangeEnd):
                #print(entry, current)
                pass
                # this never happens. Is it because the gff base values are on some sort of different scale?
                # out.write(line)

gffFh.close()
outFh.close()
