from Bio.SeqRecord import SeqRecord
import requests
from Bio import SeqIO


def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq=nuc_record.seq.translate(cds=True), \
                     id="trans_" + nuc_record.id, \
                     description="translation of CDS, using default table")


inputFile = "snpGenes.fasta"
outputFile = inputFile.replace('.fasta', '1_p.fasta')

# proteins = (make_protein_record(nuc_rec) for nuc_rec in SeqIO.parse(inputFile, "fasta"))

proteins2 = []
for nuc_rec in SeqIO.parse(inputFile, "fasta"):
    try:  # this will remove sequences that do not begin with a start codon.
        proteins2.append(make_protein_record(nuc_rec))
    except:
        continue

SeqIO.write(proteins2, outputFile, "fasta")

f_open = open(outputFile, "rU")

n = 0
for rec in SeqIO.parse(f_open, "fasta"):
    id = rec.id
    seq = rec.seq

    fasta_seq = ">" + str(id) + "\n" + str(seq)

    run = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5/run/'
    email = 'cslaubau@purdue.edu'
    title = 'projectResults'
    sequence = fasta_seq

    # send the initial query
    command = {'email': email, 'title': title, 'sequence': sequence}
    response = requests.post(run, command)
    id = response.text
    print('job {} submitted'.format(id))

    # poll for job completion
    import time

    maxtries = 1000
    notready = 1

    status = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/'

    while notready:
        response = requests.get(status + id)
        print('    polling... response->{}'.format(response.text))

        if 'FINISHED' in response.text:
            notready = 0
            break
        else:
            notready += 1

        if notready >= maxtries:
            break

        # don't poll too often
        time.sleep(20)

    if notready > 0:
        # polling reached limit
        print('unable to find result {} in {} tries'.format(id, notready))
        # exit(1)
        continue
    else:
        print('interproscan {} finished'.format(id))

        # get the final result
    result = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/'
    result += id + '/{}'.format('xml')
    response = requests.get(result)

    name = "output" + str(n) + ".xml"
    n += 1
    interpro_result = open(name, "w")
    interpro_result.write('{}\n'.format(response.text))

    # print('\n', response.text)

