"""
condenseGWAS.py takes in a GWAS 'txt' format file with columns:
    {SNP name, chromosome, position, log(p) value}

Checks each SNP log(p) value, saving those above a given threshold
Exports key snips to '.cgwas' file

Created for BIOL595 Project
"""

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
