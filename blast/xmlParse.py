"""Takes in a .xml file search result from Blast and saves the pertinent results. ####NEEDS WORK

    (script follows blastPoll.py in order)
"""

from Bio.Blast import NCBIXML as xml

geneList = 'Netblotch.genes'
fhGenes = open(geneList, 'r')
genes = []
for line in fhGenes:
    line.strip()
    m = line.split('\t')
    m[1] = m[1][0:-1]
    genes.append(m)
fhGenes.close()


input = ["snpGenes0.xml", "snpGenes1.xml", "snpGenes2.xml"]
results = []
for file in input:
    result_handle = open(file)
    blast_records = xml.parse(result_handle)
    results.append(blast_records)

# blast_records = list(blast_records)

fh = open(input[0].replace('0.xml', 'All.out'), 'w')

numSave = 5            # the number of entries to save for each search parameter
eThreshold = 10e-20     # the e-value threshold to include

entry = 0
for blast_records in results:
    for blast_record in blast_records:
        n = 0
        print(genes[entry])
        entry += 1
        print(blast_record.descriptions)
        try:
            fh.write(genes[entry][0] + '\t' + genes[entry][1] + '\n')
        except:
            break
        for desc in blast_record.descriptions:
            n += 1

            if n > numSave: # break out of the loop if the number of entries to be saved has been reached
                break

            print( desc.title)
            print('e-value =' + str(desc.e))
            print('score: ' + str(desc.score))
            #if (desc.e <= eThreshold):
                #print('^ entry would be saved ^')
                # add a way of marking which sequence this correlates with

            fh.write('\t' + str(desc) + '\n')
            # print(desc.num_alignments)

fh.close()