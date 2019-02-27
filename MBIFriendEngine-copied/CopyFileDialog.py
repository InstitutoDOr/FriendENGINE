# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'CopyFileDialog.ui'
#
# Created: Mon Aug 31 17:49:33 2015
#      by: PyQt4 UI code generator 4.10.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
import os

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def __init__(self,copyFunc):
        self.copyFunc = copyFunc
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(431, 278)
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(30, 240, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.widget = QtGui.QWidget(Dialog)
        self.widget.setGeometry(QtCore.QRect(0, 0, 431, 231))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.widget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(self.widget)
        self.label.setWordWrap(True)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.lineEdit = QtGui.QLineEdit(self.widget)
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.gridLayout.addWidget(self.lineEdit, 0, 1, 1, 1)
        self.folderBut = QtGui.QPushButton(self.widget)
        self.folderBut.setObjectName(_fromUtf8("folderBut"))
        self.gridLayout.addWidget(self.folderBut, 0, 2, 1, 1)
        self.copyBut = QtGui.QPushButton(self.widget)
        self.copyBut.setObjectName(_fromUtf8("copyBut"))
        self.gridLayout.addWidget(self.copyBut, 1, 1, 1, 1)
        self.label_2 = QtGui.QLabel(self.widget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.copyBut.setEnabled(False)
        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QObject.connect(self.folderBut, QtCore.SIGNAL(_fromUtf8("clicked()")), self.getFolderName)
        QtCore.QObject.connect(self.copyBut, QtCore.SIGNAL(_fromUtf8("clicked()")), self.copyData)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def getFolderName(self):
        print "GetFolderName"
        fname = QtGui.QFileDialog.getExistingDirectory(caption="Save directory")
        self.lineEdit.setText(fname)
        self.copyBut.setEnabled(True)
    
    def copyData(self):
        if os.access(self.lineEdit.text(),os.F_OK):
            self.copyFunc(self.lineEdit.text())

    
    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label.setText(_translate("Dialog", "Please note that this will delete all files in the NIFTI and DICOM for this study. Use the buttons below to copy data before proceeding.", None))
        self.folderBut.setText(_translate("Dialog", "Folder", None))
        self.copyBut.setText(_translate("Dialog", "Copy", None))
        self.label_2.setText(_translate("Dialog", "Local Location", None))

