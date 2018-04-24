"""
Uses the .cgwas file created in 'condenseGWAS.py' and the sequence gff file (in .txt format) as input.
Cross-references the significant GWAS SNPs in .cgwas with the .gff file.
Outputs a .genes file containing the gene detected and SNPs associated with it

(script follows csv2txt.py and condenseGWAS.py in order)
"""

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
                snpGeneOverlaps[snp].append(ID)
                #print(entry, current)
                # out.write(line)

outFile = gwasFile.replace('.cgwas', '.genes')
try:
    outFh = open(outFile, 'w')
except IOError:
    print('There was an error opening the output file')
    exit(1)

for k in snpKeys:
    for gene in snpGeneOverlaps[k]:
        outFh.write(gene + '\t' + k + '\n')

gffFh.close()
outFh.close()
