gffFile = 'Hv_IBSC_PGSB_v2p37_genes.csv'
try:
    gffFh = open(gffFile, 'r')
except IOError:
    print('There was an error opening the gff file')
    exit(1)

gffModFile = gffFile.replace('.csv', '.txt')
try:
    modFh = open(gffModFile, 'w')
except IOError:
    print('There was an error opening the output file')
    exit(1)

for line in gffFh:
    if 'gene_id' in line:
        current = line.split(',')  # need 0, 3, 4, and 8 --> [chromosome, begin, end, transcript ID]
        del current[5:9], current[1:3]

        if len(current) > 4: # trim the excess from oddly formatted lines
            while len(current) > 4:
                del current[-1]

        chrom = current[0]
        if 'Pt' in chrom:
            pass
        elif 'Un' in chrom:
            chrom = 'Un'
        else:
            chrom = chrom[3:-1]

        current[0] = chrom

        modFh.write('{}\t{}\t{}\t{}\n'.format(current[0], current[1], current[2], current[3]))

gffFh.close()
modFh.close()
