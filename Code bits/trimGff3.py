gffFile = 'Triticum_aestivum.TGACv1.38.gff3'
try:
    gffFh = open(gffFile, 'r')
except IOError:
    print('There was an error opening the gff file')
    exit(1)

gffModFile = gffFile.replace('.gff3', '.txt')
try:
    modFh = open(gffModFile, 'w')
except IOError:
    print('There was an error opening the output file')
    exit(1)

for line in gffFh:
    if line.startswith('#'):
        continue
    elif 'ID=transcript:' in line:
        current = line.split()  # need 0, 3, 4, and 8 --> [chromosome, begin, end, transcript ID]
        del current[5:8], current[1:3]
        chrom = current[0].split('_')[-1]

        transcriptID = current[3].split(':')[1]
        transcriptID = transcriptID.split(';')[0]
        transcriptID = transcriptID.split('.')[0]

        # if chrom[0].isdigit():  # trim the excess off of the numerically named chromosomes
        #     chrom = chrom[0:2]

        if len(current) > 4:  # trim the excess from oddly formatted lines
            while len(current) > 4:
                del current[-1]

        current[0] = chrom
        current[3] = transcriptID

        modFh.write('{}\t{}\t{}\t{}\n'.format(current[0], current[1], current[2], current[3]))

gffFh.close()
modFh.close()
