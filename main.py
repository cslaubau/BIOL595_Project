"""
Checks each SNP log(p) value, saving those above a given threshold
Exports key snips to '.cgwas' file

Created for BIOL595 Project
"""

# Part 1

# Open input file
inFile = 'Netblotch.txt'
try:
    inFh = open(inFile, 'r')
    print('File loaded')
except IOError:
    print('There was an error opening the GWAS input file')
    exit(1)

# Create output file
outFile = inFile.replace('.txt', '.cgwas')
try:
    outFh = open(outFile, 'w')
except IOError:
    print('There was an error opening the condensed GWAS output file')
    exit(1)

logPCutoff = 3.0  # What is the typical log(p) cutoff you want to check? Eventually make this a user input
addWidth = 2e5 # What is the typical log(p) cutoff you want to check? Eventually make this a user input
# Accidentally did double... rerun everything with 1e5

skipLine = 1  # Does the file have a header? Could be input to the user or something we check automatically
numSeqIn = 0
numSeqOut = 0
[idPrev, chromPrev, beginPrev, endPrev, logPPrev] = ["", "", "", "", ""]
hold = []
for line in inFh:
    if skipLine:
        skipLine = 0
        continue
    else:
        numSeqIn += 1
        if float(line.split()[3]) >= logPCutoff:
            numSeqOut += 1

            # Convert position and find beginning and end
            [idCurrent, chromCurrent, midCurrent, logPCurrent] = line.split()
            if len(idCurrent.split('_')) > 1:
                # S(Chrome)_## type
                #midCurrent = idCurrent.split('_')[1]
                beginCurrent = int(midCurrent) - addWidth
                if beginCurrent < 0:
                    beginCurrent = 0
                beginCurrent = str(int(beginCurrent))
                endCurrent = str(int(int(midCurrent) + addWidth))
            else:
                # Irregular
                #midCurrent = int(int(midCurrent) * 1e6)
                beginCurrent = int(midCurrent) - addWidth
                if beginCurrent < 0:
                    beginCurrent = 0
                beginCurrent = str(int(beginCurrent))
                endCurrent = str(int(int(midCurrent) + addWidth))

            # Check for overlap and combine
            if (chromPrev in chromCurrent) and (chromCurrent in chromPrev):
                if (int(beginCurrent) >= int(beginPrev)) and (int(beginCurrent) <= int(endPrev)):
                    beginCurrent = beginPrev
                    idCurrent = idPrev + ';' + idCurrent
                    del hold[-1]
                else:
                    pass
            else:
                pass

            lineCurrent = idCurrent + '\t' + chromCurrent + '\t' + beginCurrent + '\t' + endCurrent + '\t' + logPCurrent + '\n'
            [idPrev, chromPrev, beginPrev, endPrev, logPPrev] = [idCurrent, chromCurrent, beginCurrent, endCurrent, logPCurrent]

            hold.append(lineCurrent)

for entry in hold:
    outFh.write('{}'.format(entry))

inFh.close()
outFh.close()

print("\n{} sequences read in from '{}'.".format(numSeqIn, inFile))
print("{} sequences with log(p) value greater than {} output to '{}'.".format(numSeqOut, logPCutoff, outFile))

# Part 2
# Converts a .csv formatted gff file to a .txt file with only the pertinent information.

gffFile = 'Hv_IBSC_PGSB_v2p37_genes.csv'
try:
    gffFh = open(gffFile, 'r')
except IOError:
    print('There was an error opening the gff file')
    exit(1)

gffModFile = gffFile.replace('.csv', '.txt')
try:
    modFh = open(gffModFile, 'w')
except IOError:
    print('There was an error opening the output file')
    exit(1)

for line in gffFh:
    if 'gene_id' in line:
        current = line.split(',')  # need 0, 3, 4, and 8 --> [chromosome, begin, end, transcript ID]
        del current[5:9], current[1:3]

        if len(current) > 4: # trim the excess from oddly formatted lines
            while len(current) > 4:
                del current[-1]

        chrom = current[0]
        if 'Pt' in chrom:
            pass
        elif 'Un' in chrom:
            chrom = 'Un'
        else:
            chrom = chrom[3:-1]

        current[0] = chrom

        modFh.write('{}\t{}\t{}\t{}\n'.format(current[0], current[1], current[2], current[3]))

gffFh.close()
modFh.close()

# Part 3
# Uses the .cgwas file created in 'condenseGWAS.py' and the sequence gff file (in .txt format) as input.
# Cross-references the significant GWAS SNPs in .cgwas with the .gff file.
# Outputs a .genes file containing the gene detected and SNPs associated with it

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


# Part 4
# takes in the gene targets from the GWAS search (as .genes) and the full cds of the target organism (as .fa)
# For each target gene, as search is done in the cds, and the first matching sequence is saved to output (.fasta)


# Open input gene file
inGeneFile = 'Netblotch.genes'
try:
    genesFh = open(inGeneFile, 'r')
    print('File loaded')
except IOError:
    print('There was an error opening the gene input file')
    exit(1)

# Open input fasta file
fastaFile = 'Hordeum_vulgare.Hv_IBSC_PGSB_v2.cds.all.fa'
try:
    fastaFh = open(fastaFile, 'r')
    print('File loaded')
except IOError:
    print('There was an error opening the gene input file')
    exit(1)

# Create output file
outFile = 'snpGenes.fasta'
try:
    outFh = open(outFile, 'w')
except IOError:
    print('There was an error opening the Fasta output file')
    exit(1)

geneIds = []
for entry in genesFh:
    geneIds.append(entry.split()[0])
genesFh.close()

# for each gene ID
#  index the fasta file
#  find the sequence of interest
#  write to fasta output
startSaving = 0
initial = 0
for gene in geneIds:
    for line in fastaFh:
        if gene in line:
            startSaving = 1

        if startSaving and line.startswith('>'):
            initial += 1

        if startSaving and (initial > 1) and line.startswith('>'):
            fastaFh.seek(0)
            startSaving = 0
            initial = 0
            break

        if startSaving:
            outFh.write(line)

fastaFh.close()
outFh.close()

print("It finished")