# Open condensed GWAS file
gwasFile = 'water_absorption.cgwas'
try:
    gwasFh = open(gwasFile, 'r')
except IOError:
    print('There was an error opening the condensed GWAS input file')
    exit(1)

gffFile = 'Triticum_aestivum.TGACv1.38.gff3'
try:
    gffFh = open(gffFile, 'r')
except IOError:
    print('There was an error opening the gff file')
    exit(1)

# Save the ranges of interest to dictionary
gwasRanges = {} # {chrom: [begin, end]}
keys = []
for line in gwasFh:
    chrom = line.split()[1]
    if chrom not in keys:
        keys.append(chrom)
        gwasRanges[chrom] = []
    gwasRanges[chrom].append(line.split()[1:3])
gwasFh.close()

# Find and save all gff ids for ranges that lie within target ranges

gffFh.close()