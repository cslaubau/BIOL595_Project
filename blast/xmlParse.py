"""Takes in a .xml file search result from Blast and saves the pertinent results. ####NEEDS WORK

    (script follows blastPoll.py in order)
"""

from Bio.Blast import NCBIXML as xml

# result_handle = open("snpGenes.xml")
result_handle = open("test.xml")

blast_records = xml.parse(result_handle)

# blast_records = list(blast_records)

# fh = open('test.out', 'w')

numSave = 10            # the number of entries to save for each search parameter
eThreshold = 10e-20     # the e-value threshold to include

for blast_record in blast_records:
    n = 0
    for desc in blast_record.descriptions:
        n += 1

        # if n > numSave: # break out of the loop if the number of entries to be saved has been reached
        #     break

        # print(desc)
        print('\n' + desc.title)
        print('e-value =' + str(desc.e))
        print('score: ' + str(desc.score))
        if (desc.e <= eThreshold):
            print('^ entry would be saved ^')
            # add a way of marking which sequence this correlates with
            # fh.write(str(desc) + '\n')
        # print(desc.num_alignments)
