"""
This project was submitted May 4th, 2018 by Curtis Slaubaugh, Rupesh Gaire, and Apostolia Topaloudi for BIOL 595

The script takes three inputs from the user: a gwas file, a gff file for the organism, and a fasta file of the
organism's cds. (Here the example inputs are 'Netblotch.txt', 'Hv_IBSC_PGSB_v2p37_genes.csv', and
'Hordeum_vulgare.Hv_IBSC_PGSB_v2.cds.all.fa')

A Manhattan plot can be generated from the GWAS file and -log(p) cutoff.
From these inputs, the user is asked for a -log(p) cutoff for which to isolate SNPs of interest from the GWAS file.
With the SNPs of interest, the program will generate a list of base pair ranges (with the range size determined by the
user, default 200k). The gff file will then be used to isolate coding named coding sequences for which these SNP regions
overlap. The named coding sequences are then retreived from the input fasta file and saved to a new key fasta file.

With the key fasta the user can perform a Blastx and InterPro search. The Blastx search will save a user defined number
of results for each cds. In this version of the project, the local search feature does not function. The fasta file is
split into smaller files containing a maximum of 10 coding sequences in order to reserve blast computing power.
The user can define an output file to which the blast results are saved.

The InterPro search takes the key fasta file, parses individual sequences, converts the sequence to amino acid
sequences, and searches the remote InterPro database. The returned xml data is then parsed and key information about
each result saved with respect to its sequence of origin. Sequences that do not begin with a start codon are excluded
from the searches.

It is not recommended that the Blast and Interpro search are performed in the same session (call of gui.py). They are
both remote processes that take a long time to process.
"""

from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog

from UI_Files import GUI_UI
import sys, time, datetime
import subprocess as sub
from time import sleep
from Bio.Blast.Applications import NcbiblastxCommandline as cl

from fasta import fastaNext

from Bio.Blast import NCBIXML as xml
import os
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
import requests
from Bio import SeqIO

stimTimer = QtCore.QTimer()
guiTimer = QtCore.QTimer()

class MainWindow(QtWidgets.QMainWindow, GUI_UI.Ui_Form):
    """The main window is a gui object that the user interacts with."""

    def __init__(self):
        super().__init__()
        self.ui = GUI_UI.Ui_Form()
        self.ui.setupUi(self)

        self.gwasFileName = ""
        self.gffFileName = ""
        self.fastaFileName = ""
        self.filesLoaded = [0, 0, 0]
        self.logPCutoff = self.ui.logp_le.text()
        self.addWidth = self.ui.baseWidth_le.text()
        self.blastJobs = []
        self.interproJobs = []
        self.interproResultFiles = []

        self.connectActions()  # attach functions to the interactable objects in the GUI

        # Limit the user's ability to cause errors by disabling functions until prior requirements are met
        self.ui.gff_gb.setEnabled(False)
        self.ui.fasta_gb.setEnabled(False)
        self.ui.parameters_gb.setEnabled(False)

        self.ui.gwas_le.setEnabled(False)
        self.ui.gff_le.setEnabled(False)
        self.ui.fasta_le.setEnabled(False)

        self.ui.local_rb.setEnabled(False)
        self.ui.remote_rb.setChecked(True)

        self.ui.search_gb.setEnabled(False)
        self.ui.blast_gb.setEnabled(False)
        self.ui.interpro_gb.setEnabled(False)

        self.currentDateTime = datetime.datetime.now()

    def connectActions(self):  # buttons, etc
        # input -- load the associated file types
        self.ui.gwas_pb.clicked.connect(self.openGwasFile)
        self.ui.gff_pb.clicked.connect(self.openGffFile)
        self.ui.fasta_pb.clicked.connect(self.openFastaFile)

        # set parameters
        self.ui.snpSearch_pb.clicked.connect(self.mySnpSearchButtonClicked)  # run the SNP search

        # general search
        self.ui.blastEnable_cb.stateChanged.connect(self.blastEnableStateChanged)  # disable and enable blast
        self.ui.interproEnable_cb.stateChanged.connect(self.interproEnableStateChanged)  # disable and enable interpro

        # blast search
        self.ui.numSave_le.textChanged.connect(self.myNumSaveLineEditModified)  # detect errors in input
        self.ui.doSearch_pb.clicked.connect(self.searchMain)  # run the main search

        guiTimer.timeout.connect(self.updateGui)

    def openGwasFile(self):
        """Creates a window for loading gwas files"""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "Open GWAS file", "", "Text Files (*.txt);;All Files (*)",
                                                  options=options)
        if fileName:
            self.gwasFileName = fileName
            try:
                self.gwasFh = open(self.gwasFileName, 'r')
                self.ui.gwas_le.setText(str(self.gwasFileName))
                print('File loaded')
            except IOError:
                print('There was an error opening the GWAS input file')
                exit(1)
            self.filesLoaded[0] = 1
            self.ui.gff_gb.setEnabled(True)
        self.updateGui()

    def openGffFile(self):
        """Creates a window for loading gff files"""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "Open gff file", "", "GFF Files (*.csv);;All Files (*)",
                                                  options=options)
        if fileName:
            self.gffFileName = fileName
            try:
                self.gffFh = open(self.gffFileName, 'r')
                self.ui.gff_le.setText(str(self.gffFileName))
                print('File loaded')
            except IOError:
                print('There was an error opening the gff file')
                exit(1)
            self.filesLoaded[1] = 1
            self.ui.fasta_gb.setEnabled(True)
        self.updateGui()

    def openFastaFile(self):
        """Creates a window for loading fasta files"""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "Open fasta file", "", "FastA Files (*.fa);;All Files (*)",
                                                  options=options)
        if fileName:
            self.fastaFileName = fileName
            try:
                self.fastaFh = open(self.fastaFileName, 'r')
                self.ui.fasta_le.setText(str(self.fastaFileName))
                print('File loaded')
            except IOError:
                print('There was an error opening the gene input file')
                exit(1)
            self.filesLoaded[2] = 1
        self.updateGui()

    def is_number(self, s):
        """Checks if the string of interest (s) can be converted to a float without error"""
        try:
            float(s)
            return True
        except ValueError:
            return False

    def mySnpSearchButtonClicked(self):
        """Activated by the 'Find Significant SNPs' button. Starts the process of trimming the SNP file, sorting the gff
         file, and generating a targeted fasta file."""
        if self.is_number(self.ui.logp_le.text()) and self.is_number(self.ui.baseWidth_le.text()):
            self.logPCutoff = float(self.ui.logp_le.text())
            self.addWidth = float(self.ui.baseWidth_le.text())
            self.addWidth = int(self.addWidth)
            # do the snp search
            self.snpSearch()
            if self.ui.generateTree_cb.isChecked():  # if phylo tree cb checked
                self.generatePhylogenicTree()
            if self.ui.generateTree_cb_2.isChecked():  # if manhattan plot cb checked
                self.generateManhattanPlot()
            self.ui.snpSearch_pb.setEnabled(False)
            QtWidgets.QMessageBox.warning(self, 'Message', "SNP search completed.", QtWidgets.QMessageBox.Ok)
            self.ui.search_gb.setEnabled(True)
        else:  # if the logp value entered or the base pair cutoff are not numerical, tell the user to fix them.
            QtWidgets.QMessageBox.warning(self, 'Message',
                                          "You have entered a value that is not a number. Please fix this before continuing.",
                                          QtWidgets.QMessageBox.Ok)

    def snpSearch(self):
        # Part 1

        # Create output file for condensed gwas file
        self.cgwasFile = self.gwasFileName.replace('.txt', '.cgwas')
        try:
            self.cgwasFh = open(self.cgwasFile, 'w')
        except IOError:
            print('There was an error opening the condensed GWAS output file')
            exit(1)

        skipLine = 1  # Does the file have a header? Could be input to the user or something we check automatically
        numSeqIn = 0
        numSeqOut = 0
        [idPrev, chromPrev, beginPrev, endPrev, logPPrev] = ["", "", "", "", ""]
        hold = []
        for line in self.gwasFh:
            if skipLine:
                skipLine = 0
                continue
            else:
                numSeqIn += 1
                if float(line.split()[3]) >= self.logPCutoff:
                    numSeqOut += 1

                    # Convert position and find beginning and end
                    [idCurrent, chromCurrent, midCurrent, logPCurrent] = line.split()
                    if len(idCurrent.split('_')) > 1:
                        # S(Chrome)_## type
                        # midCurrent = idCurrent.split('_')[1]
                        beginCurrent = int(midCurrent) - self.addWidth
                        if beginCurrent < 0:
                            beginCurrent = 0
                        beginCurrent = str(int(beginCurrent))
                        endCurrent = str(int(int(midCurrent) + self.addWidth))
                    else:
                        # Irregular
                        # midCurrent = int(int(midCurrent) * 1e6)
                        beginCurrent = int(midCurrent) - self.addWidth
                        if beginCurrent < 0:
                            beginCurrent = 0
                        beginCurrent = str(int(beginCurrent))
                        endCurrent = str(int(int(midCurrent) + self.addWidth))

                    # Check for overlap and combine
                    if (chromPrev in chromCurrent) and (chromCurrent in chromPrev):
                        if (int(beginCurrent) >= int(beginPrev)) and (int(beginCurrent) <= int(endPrev)):
                            beginCurrent = beginPrev
                            idCurrent = idPrev + ';' + idCurrent
                            del hold[-1]
                        else:
                            pass
                    else:
                        pass

                    lineCurrent = idCurrent + '\t' + chromCurrent + '\t' + beginCurrent + '\t' + endCurrent + '\t' + logPCurrent + '\n'
                    [idPrev, chromPrev, beginPrev, endPrev, logPPrev] = [idCurrent, chromCurrent, beginCurrent,
                                                                         endCurrent, logPCurrent]

                    hold.append(lineCurrent)

        for entry in hold:
            self.cgwasFh.write('{}'.format(entry))

        self.gwasFh.close()
        self.cgwasFh.close()

        # Part 2
        # Converts a .csv formatted gff file to a .txt file with only the pertinent information.
        self.gffModFile = self.gffFileName.replace('.csv', '.txt')
        try:
            self.modFh = open(self.gffModFile, 'w')
        except IOError:
            print('There was an error opening the output file')
            exit(1)

        for line in self.gffFh:
            if 'gene_id' in line:
                current = line.split(',')  # need 0, 3, 4, and 8 --> [chromosome, begin, end, transcript ID]
                del current[5:9], current[1:3]

                if len(current) > 4:  # trim the excess from oddly formatted lines
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

                self.modFh.write('{}\t{}\t{}\t{}\n'.format(current[0], current[1], current[2], current[3]))

        self.gffFh.seek(0)
        self.modFh.close()

        # Part 3
        # Uses the .cgwas file created in 'condenseGWAS.py' and the sequence gff file (in .txt format) as input.
        # Cross-references the significant GWAS SNPs in .cgwas with the .gff file.
        # Outputs a .genes file containing the gene detected and SNPs associated with it

        # Open condensed GWAS file
        try:
            self.cgwasFh = open(self.cgwasFile, 'r')
        except IOError:
            print('There was an error opening the condensed GWAS input file')
            exit(1)

        try:
            self.modFh = open(self.gffModFile, 'r')
        except IOError:
            print('There was an error opening the output file')
            exit(1)

        # Save the ranges of interest to dictionary
        gwasRanges = {}  # {chrom: [begin, end]}
        keys = []
        for line in self.cgwasFh:
            chrom = line.split()[1]
            if chrom not in keys:
                keys.append(chrom)
                gwasRanges[chrom] = []
            lineSplit = line.split()
            del lineSplit[4], lineSplit[1]
            gwasRanges[chrom].append(lineSplit)
        self.cgwasFh.close()

        # Find and save all gff ids for ranges that lie within target ranges
        snpGeneOverlaps = {}
        snpKeys = []
        for line in self.modFh:
            current = line.split()
            if current[0] in keys:
                [chrom, gffBegin, gffEnd, ID] = current
                [gffBegin, gffEnd] = [int(gffBegin), int(gffEnd)]
                # beginning or end lies between the current gwas range, save it to file
                for entry in gwasRanges[chrom]:
                    [snp, rangeBegin, rangeEnd] = [entry[0], int(entry[1]), int(entry[2])]
                    if (gffBegin >= rangeBegin and gffBegin <= rangeEnd) or (
                            gffEnd >= rangeBegin and gffEnd <= rangeEnd):
                        if snp not in snpKeys:
                            snpKeys.append(snp)
                            snpGeneOverlaps[snp] = []
                        snpGeneOverlaps[snp].append(ID)

        self.genesFile = self.gwasFileName.replace('.txt', '.genes')
        try:
            self.genesFh = open(self.genesFile, 'w')
        except IOError:
            print('There was an error opening the output file')
            exit(1)

        for k in snpKeys:
            for gene in snpGeneOverlaps[k]:
                self.genesFh.write(gene + '\t' + k + '\n')

        self.gffFh.close()
        self.genesFh.close()

        # Part 4
        # takes in the gene targets from the GWAS search (as .genes) and the full cds of the target organism (as .fa)
        # For each target gene, as search is done in the cds, and the first matching sequence is saved to output (.fasta)

        # Open input gene file
        try:
            self.genesFh = open(self.genesFile, 'r')
            print('File loaded')
        except IOError:
            print('There was an error opening the gene input file')
            exit(1)

        # Create output file
        self.snpFile = 'snpGenes.fasta'  ###########
        try:
            self.snpFh = open(self.snpFile, 'w')
        except IOError:
            print('There was an error opening the Fasta output file')
            exit(1)

        geneIds = []
        for entry in self.genesFh:
            geneIds.append(entry.split()[0])
        self.genesFh.close()

        # for each gene ID
        #  index the fasta file
        #  find the sequence of interest
        #  write to fasta output
        startSaving = 0
        initial = 0
        for gene in geneIds:
            for line in self.fastaFh:
                if gene in line:
                    startSaving = 1

                if startSaving and line.startswith('>'):
                    initial += 1

                if startSaving and (initial > 1) and line.startswith('>'):
                    self.fastaFh.seek(0)
                    startSaving = 0
                    initial = 0
                    break

                if startSaving:
                    self.snpFh.write(line)

        self.fastaFh.close()
        self.snpFh.close()

        print("\n{} sequences read in from '{}'.".format(numSeqIn, self.gwasFileName))
        print("{} sequences with log(p) value greater than {} output to '{}'.".format(numSeqOut, self.logPCutoff,
                                                                                      self.cgwasFile))
        print("FastA sequences of interest were written to {}".format(self.snpFile))

    def generatePhylogenicTree(self):
        """Generates a phylogenic tree when the corresponding check box is selected."""
        clustalw_exe = r"./clustalw2"
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile="snpGenes.fasta")
        assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
        stdout, stderr = clustalw_cline()

        tree = Phylo.read("snpGenes.dnd", "newick")

        Phylo.draw_ascii(tree)

        try:  # This often doesn't work in pycharm. Will run best from the terminal
            tree.rooted = True
            Phylo.draw(tree)
        except:
            pass

    def generateManhattanPlot(self):
        gwasFile = self.gwasFileName
        columnFormat = "3,1,2,0"
        logp = self.logPCutoff
        numHeaderLines = 1
        cmd = "python3 manhattan_plot.py " + gwasFile + ' ' + columnFormat + ' ' + str(logp) + ' ' + str(numHeaderLines)

        job = sub.Popen(cmd, shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
        job.wait()
        print("Manhattan plot saved to {}".format('manhattan-plot-' + gwasFile.split('/')[-1].split('.')[0] + '.png'))

    def myNumSaveLineEditModified(self):
        """When the number of saved lines line edit box value is changed, check that the value entered is a number.
        If it is not, tell the user to fix it and disable the search function."""
        if self.is_number(self.ui.numSave_le.text()):
            self.ui.doBlastSearch_pb.setEnabled(True)
        elif self.ui.numSave_le.text() == "":
            self.ui.doBlastSearch_pb.setEnabled(False)
        else:
            self.ui.doBlastSearch_pb.setEnabled(False)
            QtWidgets.QMessageBox.warning(self, 'Message',
                                          "You have entered a value that is not a number. Please fix this before continuing.",
                                          QtWidgets.QMessageBox.Ok)

    def blastEnableStateChanged(self):
        """When the blast enable check box is checked, either enable or disable the blast search."""
        state = self.ui.blastEnable_cb.isChecked()
        if state:
            self.ui.blast_gb.setEnabled(True)
        else:
            self.ui.blast_gb.setEnabled(False)
        return

    def interproEnableStateChanged(self):
        """When the interpro enable check box is checked, either enable or disable the interpro search. This function is
        redundant to 'blastEnableStateChanged' only because PyQt5 has strict rules for the connecting actions."""
        state = self.ui.interproEnable_cb.isChecked()
        if state:
            self.ui.interpro_gb.setEnabled(True)
        else:
            self.ui.interpro_gb.setEnabled(False)
        return

    def searchMain(self):
        """When the 'Perform Search' button is clicked, this function checks whether blast or interpro are checked and
        runs the corresponding search(s)."""
        if self.ui.blastEnable_cb.isChecked() and self.ui.interproEnable_cb.isChecked():
            # blast setup
            self.splitFastaFile()
            print("file split")
            self.blastSearchSetup()
            print("blast complete")

            # interpro setup
            self.interprSearchSetup()

            # polling
            self.pollJobs()
            print("Search complete, parsing data...")

            # parse the output
            self.blastParse()
            print("blast result parsed")
            self.interproParse()
            print("interpro result parsed")
        elif self.ui.blastEnable_cb.isChecked():
            # blast setup
            self.splitFastaFile()
            print("file split")
            self.blastSearchSetup()
            print("blast complete")

            # polling
            self.pollJobs()
            print("Search complete, parsing data...")

            # parse the output
            self.blastParse()
            print("blast result parsed")
        elif self.ui.interproEnable_cb.isChecked():
            # interpro setup
            self.interproSearchSetup()

            # parse the output
            self.interproParse()
            print("interpro result parsed")
        else:
            pass

    def splitFastaFile(self):
        """This function chops up the important SNP fasta file into 10 sequence or less files so that they can be used
        in blast search."""
        inFile = self.snpFile
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
        self.blastFastaFileNames = []
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
                self.blastFastaFileNames.append(inFile.replace('.fasta', str(fileNo) + '.fasta'))
                fhCurrent = open(self.blastFastaFileNames[-1], 'w')
            fhCurrent.write(formatEntry(id, documentation, sequence))

            iteration += 1
        try:
            fhCurrent.close()
        except:
            pass

        inFh.close()

    def blastSearchSetup(self):
        """This script performs an online blast search on the input fasta file. Returns a .xml file that can be processed in
            Biopythons.
        """

        self.blastOutputFiles = []
        for file in self.blastFastaFileNames:
            if '.fsa' in self.blastFastaFileNames:
                self.blastOutputFiles.append(file.replace('.fa', '.xml'))
            else:
                self.blastOutputFiles.append(file.replace('.fasta', '.xml'))

        self.blastJobs = []
        for i in range(len(self.blastFastaFileNames)):
            blastx_cline = cl(query=self.blastFastaFileNames[i], db="nr", evalue=0.001, outfmt=5,
                              out=self.blastOutputFiles[i],
                              remote=1)
            job = sub.Popen(str(blastx_cline), shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
            self.blastJobs.append(job)

    def pollJobs(self):
        """This function polls until all jobs finish. It will give a runtime and the breakdown of which jobs are still
        processing. Heavily adapted from class example."""
        jobs = self.blastJobs + self.interproJobs
        numJobs = len(self.blastFastaFileNames) + 0

        done = 0
        delay = 60  # number of seconds to wait between polls
        n = numJobs
        runtime = 0
        while done < n:
            print('\nPolling')
            runtime += 1
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
                    print('\tRuntime: ' + str(runtime - 1) + ' minutes')
            sleep(delay)

    def blastParse(self):
        """Takes in a .xml file search result from Blast and saves the pertinent results."""

        geneList = self.genesFile
        try:
            fhGenes = open(geneList, 'r')
        except:
            print("There was an error reading the .genes file")
            exit(1)

        genes = []
        for line in fhGenes:
            line.strip()
            m = line.split('\t')
            m[1] = m[1][0:-1]
            genes.append(m)
        fhGenes.close()

        xmlInput = self.blastOutputFiles
        results = []
        for file in xmlInput:
            result_handle = open(file)
            blast_records = xml.parse(result_handle)
            results.append(blast_records)

        fh = open(self.ui.blastOutputFileName_le.text(), 'w')

        self.numSave = int(self.ui.numSave_le.text())  # the number of entries to save for each search parameter
        # eThreshold = 10e-20     # the e-value threshold to include

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

                    if n > self.numSave:  # break out of the loop if the number of entries to be saved has been reached
                        break

                    # print(desc.title)
                    # print('e-value =' + str(desc.e))
                    # print('score: ' + str(desc.score))
                    # if (desc.e <= eThreshold):
                    # print('^ entry would be saved ^')
                    # add a way of marking which sequence this correlates with

                    fh.write('\t' + str(desc) + '\n')

        fh.close()

        # delete mess files
        messFiles = self.blastOutputFiles + self.blastFastaFileNames

        ## if file exists, delete it ##
        for myfile in messFiles:
            if os.path.isfile(myfile):
                os.remove(myfile)
            else:  ## Show an error ##
                print("Error: %s file not found" % myfile)

    def interproSearchSetup(self):
        """Sets up and runs a remote interpro search. First splits the fasta file into single sequences and converts
        them to amino acid sequences. Discards any sequence that does not begin with a start codon. Heavily adapted
        from class example."""

        def make_protein_record(nuc_record):
            """Returns a new SeqRecord with the translated sequence (default table)."""
            return SeqRecord(seq=nuc_record.seq.translate(cds=True), id="trans_" + nuc_record.id,
                             description="translation of CDS, using default table")

        inputFile = self.snpFile
        outputFile = inputFile.replace('.fasta', '1_p.fasta')

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
            self.interproResultFiles.append(name)

            n += 1
            interpro_result = open(name, "w")
            interpro_result.write('{}\n'.format(response.text))

    def interproParse(self):
        """Parses the interpro result xml files. Saves pertinent matching infromation like accession ID, library, name,
        description. Outputs the data to user defined text file and deletes the xml files."""
        files = self.interproResultFiles
        outFile = self.ui.interproOutputFileName_le.text()
        outFh = open(outFile, 'w')

        for file in files:
            fh = open(file, 'r')

            lineNum = 1
            lines2save = []
            xref = ""
            for line in fh:
                if line.lstrip().startswith('<xref'):
                    xref = line
                if line.lstrip().startswith('<signature '):
                    lines2save.append([lineNum, 'start'])
                elif line.lstrip().startswith('</signature'):
                    lines2save.append([lineNum, 'stop'])

                lineNum += 1
            fh.seek(0)

            xref = xref.replace('<xref ', '')
            xref = xref.lstrip()
            xref = xref.replace('/>', '')
            xref = 'Gene: ' + xref

            allLines = []
            n = 0
            record = 0
            for num in range(lineNum):
                for value in lines2save:
                    if value[0] == num:
                        if value[1] == 'start':
                            record = 1
                        elif value[1] == 'stop':
                            record = 0
                if record:
                    allLines.append(n)

                n += 1

            string = ""
            k = 1
            for line in fh:
                if k in allLines:
                    string += line.lstrip()
                k += 1

            info = ""
            i = 1
            outFh.write(xref)
            for line in string.split('\n'):
                if line.startswith('<signature ac='):
                    info = line.replace('<signature ', '')
                    info = info.replace('>', '')
                    if i != 1:
                        info = '\n\t' + info
                elif line.startswith('<signature-library-release'):
                    info = line.replace('<signature-library-release ', '')
                    info = info.replace('/>', '')
                elif line.startswith('<entry ac='):
                    info = line.replace('<entry ', '')
                    info = info.replace('/>', '')

                if info != '':
                    outFh.write('\t' + info + '\n')
                info = ""
                i += 1

        outFh.close()

        # delete mess files
        messFiles = self.interproResultFiles

        ## if file exists, delete it ##
        for myfile in messFiles:
            if os.path.isfile(myfile):
                os.remove(myfile)
            else:  ## Show an error ##
                print("Error: %s file not found" % myfile)

    def updateGui(self):
        """"A function that is activated by transitional actions in the GUI. Enables or disables sections and features
        as needed."""
        inputComplete = 1  # once the input files have been loaded, allow the user to progress.
        for i in self.filesLoaded:
            inputComplete *= i
        if inputComplete:
            self.ui.parameters_gb.setEnabled(True)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())
