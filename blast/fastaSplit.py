from fasta import fastaNext

inFile = 'snpGenes.fasta'
inFh = open(inFile, 'r')

def formatEntry(id, doc, seq):
    lineLength = 60  # assign default length

    string = '>' + id + ' ' + doc + '\n'  # create the header line
    i = 0
    while i < len(seq):  # format the sequence to the specified line length
        string += seq[i:i + lineLength] + '\n'
        i += lineLength
    return string

buffer = 0
setSize = 10
iteration = 0
fileNames = []
while True:
    [id, documentation, sequence, buffer] = fastaNext(inFh, buffer)
    if not id:
        break
    if iteration % 10 == 0:
        fileNo = int(iteration // 10)
        print(fileNo)
        if iteration != 0:
            fhCurrent.close()
            # close the previous file
            pass
        fileNames.append(inFile.replace('.fasta', str(fileNo) + '.fasta'))
        fhCurrent = open(fileNames[-1], 'w')
    fhCurrent.write(formatEntry(id, documentation, sequence))

    iteration += 1
try:
    fhCurrent.close()
except:
    pass

print(fileNames)

inFh.close()