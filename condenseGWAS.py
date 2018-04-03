"""
condenseGWAS.py takes in a GWAS 'txt' format file with columns:
    {SNP name, chromosome, position, log(p) value}

Checks each SNP log(p) value, saving those above a given threshold
Exports key snips to '.cgwas' file

Created for BIOL595 Project
"""

# Open input file
inFile = 'water_absorption.txt'
try:
    inFh = open(inFile, 'r')
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

logPCutoff = 4.0  # What is the typical log(p) cutoff you want to check? Eventually make this a user input

skipLine = 1  # Does the file have a header? Could be input to the user or something we check automatically
numSeqIn = 0
numSeqOut = 0
for line in inFh:
    if skipLine:
        skipLine = 0
        continue
    else:
        numSeqIn += 1
        if float(line.split()[3]) >= logPCutoff:
            numSeqOut += 1

            # Convert position?
            # Should there be a check for redundancy?

            outFh.write('{}'.format(line))

print("\n{} sequences read in from '{}'.".format(numSeqIn, inFile))
print("{} sequences with log(p) value greater than {} output to '{}'.".format(numSeqOut, logPCutoff, outFile))

inFh.close()
outFh.close()
