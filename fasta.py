"""=================================================================================================
fasta.py
module to read each fasta sequence from a multiple fasta file
2 February 2018     Michael Gribskov
================================================================================================="""


def fastaNext(fh, buffer):
    """---------------------------------------------------------------------------------------------
    return the next entry from an open file into the object.  The only way to tell you have reached
    the end of a sequence is when you read the first line of the next file.  This line is returned
    in buffer.
    usage
        id, doc, seq, buffer = fastaNext(fasta_in, buffer)
           ...
   :param: fh - an open readable filehandle for a Fasta file
   :param buffer: the next line of the file
   :return: id, documentation, sequence, buffer
    ---------------------------------------------------------------------------------------------"""
    id = ''
    documentation = ''
    sequence = ''

    #  if buffer is empty read a line
    if buffer:
        line = buffer
        buffer = ''
    else:
        line = ''
        for line in fh:
            if line.isspace():
                continue
            else:
                break

    # not successful in finding a header line, must be end of file
    if not line:
        return id, documentation, sequence, buffer

    # get the ID and documentation from the doc line
    line = line.rstrip()
    try:
        id, documentation = line.split(" ", maxsplit=1)
    except ValueError:
        # if documentation is missing, split fails
        # print('fastaNext - documentation is missing')
        id = line

    id = id.lstrip('> ')

    # read the sequence, since the id and documentation are already parsed, it doesn't need to be
    # done here
    for line in fh:
        if line.isspace():
            # skip blank lines
            continue

        line = line.rstrip()  # remove newline
        # remove N and *
        line = line.replace('N', '')
        line = line.replace('*', '')

        if line.startswith('>'):
            # start of next sequence
            buffer = line
            break

        else:
            sequence += line

    return id, documentation, sequence, buffer

    # End of fastaNext
