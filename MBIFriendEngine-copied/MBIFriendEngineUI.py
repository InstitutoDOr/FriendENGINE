# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MBIFriendEngine.ui'
#
# Created: Fri Aug 28 11:42:31 2015
#      by: PyQt4 UI code generator 4.10.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

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

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1220, 585)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.resetParamButton = QtGui.QPushButton(self.centralwidget)
        self.resetParamButton.setObjectName(_fromUtf8("resetParamButton"))
        self.horizontalLayout.addWidget(self.resetParamButton)
        self.startButton = QtGui.QPushButton(self.centralwidget)
        self.startButton.setObjectName(_fromUtf8("startButton"))
        self.horizontalLayout.addWidget(self.startButton)
        self.stopButton = QtGui.QPushButton(self.centralwidget)
        self.stopButton.setObjectName(_fromUtf8("stopButton"))
        self.horizontalLayout.addWidget(self.stopButton)
        self.clearRunButton = QtGui.QPushButton(self.centralwidget)
        self.clearRunButton.setObjectName(_fromUtf8("clearRunButton"))
        self.horizontalLayout.addWidget(self.clearRunButton)
        self.clearButton = QtGui.QPushButton(self.centralwidget)
        self.clearButton.setObjectName(_fromUtf8("clearButton"))
        self.horizontalLayout.addWidget(self.clearButton)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.verticalLayout.addWidget(self.line)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setSpacing(-1)
        self.horizontalLayout_4.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.horizontalLayout_4.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.anatIndicator = QLed(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.anatIndicator.sizePolicy().hasHeightForWidth())
        self.anatIndicator.setSizePolicy(sizePolicy)
        self.anatIndicator.setObjectName(_fromUtf8("anatIndicator"))
        self.horizontalLayout_4.addWidget(self.anatIndicator)
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout_4.addWidget(self.label_2)
        self.ledIndicator_3 = QLed(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ledIndicator_3.sizePolicy().hasHeightForWidth())
        self.ledIndicator_3.setSizePolicy(sizePolicy)
        self.ledIndicator_3.setObjectName(_fromUtf8("ledIndicator_3"))
        self.horizontalLayout_4.addWidget(self.ledIndicator_3)
        self.funcRefIndicator = QtGui.QLabel(self.centralwidget)
        self.funcRefIndicator.setObjectName(_fromUtf8("funcRefIndicator"))
        self.horizontalLayout_4.addWidget(self.funcRefIndicator)
        self.ledIndicator = QLed(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ledIndicator.sizePolicy().hasHeightForWidth())
        self.ledIndicator.setSizePolicy(sizePolicy)
        self.ledIndicator.setObjectName(_fromUtf8("ledIndicator"))
        self.horizontalLayout_4.addWidget(self.ledIndicator)
        self.funcRunIndicator = QtGui.QLabel(self.centralwidget)
        self.funcRunIndicator.setObjectName(_fromUtf8("funcRunIndicator"))
        self.horizontalLayout_4.addWidget(self.funcRunIndicator)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.line_3 = QtGui.QFrame(self.centralwidget)
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.verticalLayout.addWidget(self.line_3)
        self.engineOutput = QtGui.QTextEdit(self.centralwidget)
        self.engineOutput.setReadOnly(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(6)
        sizePolicy.setHeightForWidth(self.engineOutput.sizePolicy().hasHeightForWidth())
        self.engineOutput.setSizePolicy(sizePolicy)
        self.engineOutput.setObjectName(_fromUtf8("engineOutput"))
        self.verticalLayout.addWidget(self.engineOutput)
        self.convOutput = QtGui.QTextEdit(self.centralwidget)
        self.convOutput.setReadOnly(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.convOutput.sizePolicy().hasHeightForWidth())
        self.convOutput.setSizePolicy(sizePolicy)
        self.convOutput.setObjectName(_fromUtf8("convOutput"))
        self.verticalLayout.addWidget(self.convOutput)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1220, 22))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MBIFriendEngine", None))
        self.startButton.setText(_translate("MainWindow", "Start", None))
        self.resetParamButton.setText(_translate("MainWindow", "Reset Params", None))
        self.stopButton.setText(_translate("MainWindow", "Stop", None))
        self.clearRunButton.setText(_translate("MainWindow", "Clear Run", None))
        self.clearButton.setText(_translate("MainWindow", "Clear", None))
        self.anatIndicator.setToolTip(_translate("MainWindow", "Led indicator/button", None))
        self.anatIndicator.setWhatsThis(_translate("MainWindow", "Led indicator/button", None))
        self.label_2.setText(_translate("MainWindow", "Anatomical Image", None))
        self.ledIndicator_3.setToolTip(_translate("MainWindow", "Led indicator/button", None))
        self.ledIndicator_3.setWhatsThis(_translate("MainWindow", "Led indicator/button", None))
        self.funcRefIndicator.setText(_translate("MainWindow", "Functional Reference Image", None))
        self.ledIndicator.setToolTip(_translate("MainWindow", "Led indicator/button", None))
        self.ledIndicator.setWhatsThis(_translate("MainWindow", "Led indicator/button", None))
        self.funcRunIndicator.setText(_translate("MainWindow", "Functional Current Run", None))

from QLed import QLed
