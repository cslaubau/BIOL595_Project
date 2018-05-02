import PyQt5

def selectFile():
    lineEdit.setText(QFileDialog.getOpenFileName())

pushButton.clicked.connect(selectFile)