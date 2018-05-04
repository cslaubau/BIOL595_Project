"""=================================================================================================
fasta.py
write procedural code to read each fasta sequence from short.fa.
report the ID and length of each sequence.

2 February 2018     Michael Gribskov
================================================================================================="""
from blast.fasta import fastaNext

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
infile = 'Triticum_aestivum.TGACv1.dna.toplevel.fa'
fasta_in = open(infile, 'r')

outfile = 'task3.' + infile + '.out'
fasta_out = open(outfile, 'w')

buffer = ''
while True:
    id, doc, seq, buffer, doWrite = fastaNext(fasta_in, buffer)
    if not id:
        break

    if doWrite:
        fasta_out.write('{} {}\n'.format(id, len(seq)))

exit(0)
