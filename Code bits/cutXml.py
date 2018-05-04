files = ["output0.xml", "output1.xml", "output2.xml", "output3.xml", "output4.xml"]
outFile = 'test.txt'
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

    # print(allLines)

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