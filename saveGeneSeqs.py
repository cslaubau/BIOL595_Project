"""
snpGeneSeqs.py takes in the gene targets from the GWAS search (as .genes) and the full cds of the target organism (as .fa)
For each target gene, as search is done in the cds, and the first matching sequence is saved to output (.fasta)

(script follows findGffEntries.py in order)
"""

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