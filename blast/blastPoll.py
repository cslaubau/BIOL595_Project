"""This script performs an online blast search on the input fasta file. Returns a .xml file that can be processed in
    Biopython (see xmlParse.py).

    (script follows saveGeneSeqs.py in order)
"""

import subprocess as sub
from time import sleep
from Bio.Blast.Applications import NcbiblastxCommandline as cl

fastaFile = 'snpGenes.fasta'
# fastaFile = 'test.fasta'

if '.fsa' in fastaFile:
    outputFile = fastaFile.replace('.fsa', '.xml')
else:
    outputFile = fastaFile.replace('.fasta', '.xml')

blastx_cline = cl(query=fastaFile, db="nr", evalue=0.001, outfmt=5, out=outputFile, remote=1)

jobs = [sub.Popen(str(blastx_cline), shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)]


# poll until all jobs finish
done = 0
delay = 60  # number of seconds to wait between polls
n = 1
runtime = 0
while done < n:
    print('\nPolling')
    for i in range(n):
        if jobs[i] == 'Done':
            continue

        print('    job {} ...'.format(i), end='')
        result = jobs[i].poll()

        if result != None:
            print('finished')
            jobs[i] = 'Done'
            done += 1

        else:
            print('still running')
            runtime += 1
            print('\tRuntime: ' + str(runtime-1) + ' minutes')

    sleep(delay)