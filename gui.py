from PyQt5 import QtGui, QtCore, QtWidgets

from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog
from PyQt5.QtGui import QIcon

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

stimTimer = QtCore.QTimer()
guiTimer = QtCore.QTimer()


class MainWindow(QtWidgets.QMainWindow, GUI_UI.Ui_Form):
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

        self.connectActions()

        self.ui.gff_gb.setEnabled(False)
        self.ui.fasta_gb.setEnabled(False)
        self.ui.parameters_gb.setEnabled(False)

        self.ui.gwas_le.setEnabled(False)
        self.ui.gff_le.setEnabled(False)
        self.ui.fasta_le.setEnabled(False)

        self.ui.local_rb.setEnabled(False)
        self.ui.remote_rb.setChecked(True)

        # Setup Trial Information
        self.currentDateTime = datetime.datetime.now()

        # self.ui.currentDate_de.setDate(self.currentDateTime)
        # self.ui.currentDate_de.setDisplayFormat('yyyy/MM/dd')
        # self.ui.currentTime_te.setTime(self.currentDateTime.time())
        # self.ui.currentTime_te.setDisplayFormat('hh:mm:ss')

        # self.ui.gb.setEnabled(False)
        # self.ui.le.setEnabled(False)
        # self.ui.cb.setEnabled(False)
        # self.ui.pb.setEnabled(False)

    def connectActions(self):  # buttons, etc
        # self.ui.pb.clicked.connect(self.myButtonClicked)
        # self.ui.le.textChanged.connect(self.myLineEditModified)
        # self.ui.cb.currentIndexChanged.connect(self.myIndexChanged)

        # input
        self.ui.gwas_pb.clicked.connect(self.openGwasFile)
        self.ui.gff_pb.clicked.connect(self.openGffFile)
        self.ui.fasta_pb.clicked.connect(self.openFastaFile)

        # set parameters
        self.ui.snpSearch_pb.clicked.connect(self.mySnpSearchButtonClicked)  # run the SNP search

        # blast search
        self.ui.numSave_le.textChanged.connect(self.myNumSaveLineEditModified)
        self.ui.doBlastSearch_pb.clicked.connect(self.blastSearchMain)

        guiTimer.timeout.connect(self.updateGui)

    def openGwasFile(self):
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
        try:
            float(s)
            return True
        except ValueError:
            return False

    def mySnpSearchButtonClicked(self):
        if self.is_number(self.ui.logp_le.text()) and self.is_number(self.ui.baseWidth_le.text()):
            print("values are numbers")
            self.logPCutoff = float(self.ui.logp_le.text())
            self.addWidth = float(self.ui.baseWidth_le.text())
            self.addWidth = int(self.addWidth)
            # do the snp search
            self.snpSearch()
            if self.ui.generateTree_cb.isChecked():
                self.generatePhylogenicTree()
            self.ui.snpSearch_pb.setEnabled(False)
        else:
            print("number error")

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
                    if (gffBegin >= rangeBegin and gffBegin <= rangeEnd) or (gffEnd >= rangeBegin and gffEnd <= rangeEnd):
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
        clustalw_exe = r"./clustalw2"
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile="snpGenes.fasta")
        assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
        stdout, stderr = clustalw_cline()

        tree = Phylo.read("snpGenes.dnd", "newick")

        Phylo.draw_ascii(tree)
        # tree.rooted = True
        # Phylo.draw(tree)

    def myNumSaveLineEditModified(self):
        if self.is_number(self.ui.numSave_le.text()):
            self.ui.doBlastSearch_pb.setEnabled(True)
        elif self.ui.numSave_le.text() == "":
            self.ui.doBlastSearch_pb.setEnabled(False)
        else:
            self.ui.doBlastSearch_pb.setEnabled(False)
            QtWidgets.QMessageBox.warning(self, 'Message',
                                          "You have entered a value that is not a number. Please fix this before continuing.",
                                          QtWidgets.QMessageBox.Ok)

    def blastSearchMain(self):
        self.splitFastaFile()
        #self.blastSearch()
        #self.blastParse()

    def splitFastaFile(self):
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
                fhCurrent = open(fileNames[-1], 'w')
            fhCurrent.write(formatEntry(id, documentation, sequence))

            iteration += 1
        try:
            fhCurrent.close()
        except:
            pass

        inFh.close()

    def blastSearch(self):
        """This script performs an online blast search on the input fasta file. Returns a .xml file that can be processed in
            Biopython (see xmlParse.py).

            (script follows saveGeneSeqs.py in order)
        """

        fastaFile = self.blastFastaFileNames

        self.blastOutputFiles = []
        for file in fastaFile:
            if '.fsa' in fastaFile:
                self.blastOutputFiles.append(file.replace('.fa', '.xml'))
            else:
                self.blastOutputFiles.append(file.replace('.fasta', '.xml'))

        # blastx_cline = cl(query=fastaFile, db="nr", evalue=0.001, outfmt=5, out=outputFile, remote=1)
        # jobs = [sub.Popen(str(blastx_cline), shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)]

        jobs = []
        for i in range(len(fastaFile)):
            blastx_cline = cl(query=fastaFile[i], db="nr", evalue=0.001, outfmt=5, out=self.blastOutputFiles[i],
                              remote=1)
            job = sub.Popen(str(blastx_cline), shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
            jobs.append(job)

        # poll until all jobs finish
        done = 0
        delay = 60  # number of seconds to wait between polls
        n = len(fastaFile)
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
                    print('\tRuntime: ' + str(runtime - 1) + ' minutes')

            sleep(delay)

    def blastParse(self):
        """Takes in a .xml file search result from Blast and saves the pertinent results. ####NEEDS WORK

            (script follows blastPoll.py in order)
        """

        geneList = 'Netblotch.genes'  ###
        fhGenes = open(geneList, 'r')
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

        # blast_records = list(blast_records)

        fh = open(xmlInput[0].replace('0.xml', 'All.out'), 'w')

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

                    print(desc.title)
                    print('e-value =' + str(desc.e))
                    print('score: ' + str(desc.score))
                    # if (desc.e <= eThreshold):
                    # print('^ entry would be saved ^')
                    # add a way of marking which sequence this correlates with

                    fh.write('\t' + str(desc) + '\n')
                    # print(desc.num_alignments)

        fh.close()

    def myButtonClicked(self):
        pass

    def myLineEditModified(self):
        pass

    def myIndexChanged(self):
        pass

    def messageBox(self):
        strMsg = "Would you like to do a thing?"
        buttonReply = QtWidgets.QMessageBox.question(self, 'Message', strMsg,
                                                     QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                     QtWidgets.QMessageBox.No)
        if buttonReply == QtWidgets.QMessageBox.Yes:
            print('Yes clicked.')
        else:
            print('No clicked.')

    def updateGui(self):
        inputComplete = 1
        for i in self.filesLoaded:
            inputComplete *= i
        if inputComplete:
            self.ui.parameters_gb.setEnabled(True)

    def recordToTable(self):
        pass
    # for entry in self.currentEntry:
    # 	currentRowCount = self.ui.Iop_tableWidget.rowCount()
    # 	currentRowCount += 1

    # 	self.ui.Iop_tableWidget.setRowCount(currentRowCount)
    # 	item = QtWidgets.QTableWidgetItem()

    # 	self.ui.Iop_tableWidget.setVerticalHeaderItem(currentRowCount-1, item)
    # 	item = self.ui.Iop_tableWidget.verticalHeaderItem(currentRowCount-1)
    # 	item.setText(str(currentRowCount))

    # 	[prePost, typeOfMeasure, eyeMeasured, iopReturned, timeReturned] = entry

    # 	item = QtWidgets.QTableWidgetItem()
    # 	item.setTextAlignment(QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
    # 	self.ui.Iop_tableWidget.setItem(currentRowCount-1, 0, item)
    # 	self.ui.Iop_tableWidget.item(currentRowCount-1,0).setText(iopReturned + ' mmHg')

    # 	item = QtWidgets.QTableWidgetItem()
    # 	item.setTextAlignment(QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
    # 	self.ui.Iop_tableWidget.setItem(currentRowCount-1, 1, item)
    # 	self.ui.Iop_tableWidget.item(currentRowCount-1,1).setText(timeReturned)

    # 	item = QtWidgets.QTableWidgetItem()
    # 	item.setTextAlignment(QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
    # 	self.ui.Iop_tableWidget.setItem(currentRowCount-1, 2, item)
    # 	self.ui.Iop_tableWidget.item(currentRowCount-1,2).setText(prePost)

    # 	item = QtWidgets.QTableWidgetItem()
    # 	item.setTextAlignment(QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
    # 	self.ui.Iop_tableWidget.setItem(currentRowCount-1, 3, item)
    # 	self.ui.Iop_tableWidget.item(currentRowCount-1,3).setText(typeOfMeasure)

    # 	item = QtWidgets.QTableWidgetItem()
    # 	item.setTextAlignment(QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
    # 	self.ui.Iop_tableWidget.setItem(currentRowCount-1, 4, item)
    # 	self.ui.Iop_tableWidget.item(currentRowCount-1,4).setText(eyeMeasured)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    # app.setWindowIcon(QtGui.QIcon('splash.png'))
    # splash_pix = QtGui.QPixmap('splash.png')
    # splash = QtWidgets.QSplashScreen(splash_pix, QtCore.Qt.WindowStaysOnTopHint)
    # splash.setMask(splash_pix.mask())
    # splash.show()
    # app.processEvents()
    # time.sleep(2)
    # del(splash)
    ex = MainWindow()
    ex.show()

    # winPreIop = EnterPreIopWindow()

    sys.exit(app.exec_())
