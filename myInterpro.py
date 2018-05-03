import subprocess as sub
from time import sleep

# sequence = """>trans_HORVU4Hr1G004820 translation of CDS, using default table
# MARLAPKAKVLRDGRWTEEEAAVLVPGDIISIKLGDIIPADARLLDGDPLRIDQSALTGE
# SLPATKGPGDGVYSGSTVKQGEIEAVVIATGVHTFFGKAAHLVDSTNQVGHFQQVLTAIG
# NFCICSIAVGMFIEIIVMYPIQHRAYRPGIDNLLVLLIGGIPIAMPTVLSVTMAIGSHRL
# SQQGAITKRMTAIEEMAGMDVLCSDKTGTLTLNKLSVDKNLVEVFEKGVTQDQVILMAAR
# ASRIENQDAIDTAIVGMLGDPKEARAGIQEVHFLPFNPTDKRTALTYIDGDGKMYRVSKG
# APEQILNLAYNKSEIAQKVHTVIDKFAERGLRSLGVAYQDVPDGRKESPGSPWHFVALLP
# LFDPPRHDSAETIERALNLGVNVKMITGDQLAIGKETGRRLGMGTNMYPSSALLGQNKDE
# SIADLPVDDLIEKADGFAGVFPEHKYEIVKRLQARKHICGMTGDGVNDAPALKKADIGIA
# VADATDAARSASDIVLTEPGLSVIISAVLTSRAIFQRMTCLQIYAVSITIRIVLGFMLLA
# LIWEFDFPPFMVLIIAILNDGTIMTISKDRVKPSPLPDSWKLAEIFTTGVVLGGYLAMMT
# VIFFWAAYKTNFFPRVFHVRSLEKTAQDDFNKMLASAVYLQVSTISQALIFVTRSRSWSF
# LERPGFLLVFAFFVAQLIATLIAVYADWAFTSIKGIGWGWAGIVWLYNLVFYFPLDIIKF
# FIRYALSGKAWDLVINQRIAFTRKKHFGKEERELKWAHAQRTLHGLQPPDAKLFPEKAGY
# NELNQMAEEAKRRAEIARLRELHTLKGHVESVVKLKGLDIDTIQQSYTV"""
#
# title = sequence.split('\n')[0]
# sequence = sequence[len(title):len(sequence)-1]
# sequence.replace('\n', '')
# print(sequence)

sequence = "testProt.fasta"
email = "cslaubau@purdue.edu"
title = "test"
outFileName = "test.gff3"
format = "GFF3"
url = "http://www.ebi.ac.uk/Tools/services/rest/iprscan5/"

baseCmd = "python3 iprscan5_urllib3.py "
baseCmd += "--email=" + email + " "
baseCmd += "--title=" + title + " "
baseCmd += "--outfile=" + outFileName + " "
baseCmd += "--outformat=" + format + " "
baseCmd += "--baseURL=" + url + " "
baseCmd += "--sequence=" + sequence

# --sequence=snpGenes.fasta
# --email=cslaubau@purdue.edu
# --title=test
# --outfile=test.gff3
# --outformat=GFF3
# --baseURL=http://www.ebi.ac.uk/Tools/services/rest/iprscan5/run/

print(baseCmd)

jobs = []
numSeqs = 1
for i in range(numSeqs):
    job = sub.Popen(str(baseCmd), shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
    jobs.append(job)

# poll until all jobs finish
done = 0
delay = 60  # number of seconds to wait between polls
n = numSeqs
runtime = 0
while done < n:
    print('\nPolling')
    for i in range(n):
        runtime += 1
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
            print('\tRuntime: ' + str(runtime-1) + ' minutes')

    sleep(delay)