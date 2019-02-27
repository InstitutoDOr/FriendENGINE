import sys
from PyQt4.QtGui import QApplication, QMainWindow, QDialog, QMessageBox, QPushButton
from MBIFriendEngineUI import Ui_MainWindow  # here you need to correct the names
from RemoteRunner import RemoteRunner
from PyQt4 import QtCore
import threading
import time
import os
import CopyFileDialog
from FileCopyProgress import FileCopyProgress

from paramiko import SSHClient
from scp import SCPClient

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

ui = ""

root = "X:\\S1\\"

ssh_key = "C:\\key\\id_rsa_nfb"

friend_server = None
converter = None
app = None
run_checks = True

def runChecks():

    global run_checks
    while run_checks:
        if os.access(os.path.join(root,"RAI.nii"),os.F_OK):
            toggleRAI(True)
        else:
            toggleRAI(False)
            
        if os.access(os.path.join(root,"RFI.nii"),os.F_OK):
            toggleRFI(True)
        else:
            toggleRFI(False)

	if os.access(os.path.join(root,"CUR_RUN"),os.F_OK):            
            if len(os.listdir(os.path.join(root,"CUR_RUN"))) > 0:
                toggleCurRun(True)
            else:
                toggleCurRun(False)

        time.sleep(1)

@QtCore.pyqtSlot()
def appendEngText(text, err=False):
 
    if err:   
        ui.engineOutput.setTextColor(QtCore.Qt.red)
    else:
        ui.engineOutput.setTextColor(QtCore.Qt.black)
    ui.engineOutput.append(text)

@QtCore.pyqtSlot()
def appendConvText(text,err=False):
    
    if err:   
        ui.convOutput.setTextColor(QtCore.Qt.red)
    else:
        ui.convOutput.setTextColor(QtCore.Qt.black)
    ui.convOutput.append(text)

def handleReset():
#    reply = QMessageBox.question( None,'Message',"Are you sure you want to reset the parameter file?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
    msgBox = QMessageBox()
    msgBox.setText("Are you sure you want to reset the parameter file?")
    msgBox.addButton(_fromUtf8('Yes'), QMessageBox.YesRole)
    msgBox.addButton(_fromUtf8('Upload'), QMessageBox.YesRole)
    msgBox.addButton(_fromUtf8('No'), QMessageBox.RejectRole)
    reply = msgBox.exec_();

    print reply
    if reply == 0:
        return
        copy = RemoteRunner("cp /usr/local/FriendENGINE/Friend_Engine_Sources/Friend_Engine_Sources/Application/study_params.txt.temp /usr/local/FriendENGINE/Friend_Engine_Sources/Friend_Engine_Sources/Application/study_params.txt",[""],ssh_key)   
        copy.sig_handler.updateInfoSig.connect(appendEngText)
        copy.start()
        time.sleep(2)   
        copy.stopExe()
    elif reply == 1:
        ## 
        d_ui.setupUi(dialog)
        result = dialog.exec_() 
        print result 

def handleStart():
    ## make 
    ui.startButton.setEnabled(False)
    ui.resetParamButton.setEnabled(False)
    ui.stopButton.setEnabled(True)

    global converter
    global friend_server
    ## start the friendserver and the converter
    friend_server = RemoteRunner("/usr/local/FriendENGINE/Friend_Engine_Sources/Friend_Engine_Sources/Application/engine.sh",[""],ssh_key)
    friend_server.sig_handler.updateInfoSig.connect(appendEngText)   
    friend_server.start() 
    converter = RemoteRunner("python",["/usr/local/SystemScripts/NFBScripts/dicom_converter.py"],ssh_key)   
    converter.sig_handler.updateInfoSig.connect(appendConvText)   
    converter.start() 
    ## pipe the outputs to the browser windows
    

def handleStop():
        
    ## make
     
    ui.startButton.setEnabled(True)
    ui.stopButton.setEnabled(False)
    ui.clearButton.setEnabled(True)

    global converter
    global friend_server
    converter.stopExe()
    friend_server.stopExe()

def handleClearRun():

    for dd in os.listdir(root+"CUR_RUN"):
        print dd
        os.remove(os.path.join(root,"CUR_RUN",dd))
    ## just delete the symlinks from the CUR_RUN folder

def copyParamFile(localFile):

    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    #client.load_host_keys(os.path.expanduser('~/.ssh/known_hosts'))
    k = paramiko.RSAKey.from_private_key_file(self.key)
    client.connect("carbon2.mbi.monash.edu", username="nfb", pkey=ssh_key)
    scp = SCPClient(client.get_transport())
    scp.put(result,"/usr/local/FriendENGINE/Friend_Engine_Sources/Friend_Engine_Sources/Application/study_params.txt")
    scp.close()

def copyData(dest):
    global app
    FileCopyProgress(app=app,src=root,dest=dest)
    ##
    
def handleClear():
    ## need a dialog to copy the data

    # pop up a dialog box 
    dialog = QDialog()
    d_ui = CopyFileDialog.Ui_Dialog(copyData)
    d_ui.setupUi(dialog)
    result = dialog.exec_() 
    ## delete stuff
    ui.resetParamButton.setEnabled(True)


def toggleRAI(on):
    ui.anatIndicator.setValue(on)
def toggleRFI(on):
    ui.ledIndicator_3.setValue(on)
def toggleCurRun(on):
    ui.ledIndicator.setValue(on)

if __name__ == "__main__":

    global app
    app = QApplication(sys.argv)
    window = QMainWindow()
    ui = Ui_MainWindow()

    ui.setupUi(window)
    ## blank everything out
    ui.resetParamButton.setEnabled(True)
    ui.startButton.setEnabled(True)
    ui.stopButton.setEnabled(False)
    ui.clearRunButton.setEnabled(True)
    ui.clearButton.setEnabled(False)
    
    ## connect up the buttons
    ui.resetParamButton.clicked.connect(handleReset)
    ui.startButton.clicked.connect(handleStart)
    ui.stopButton.clicked.connect(handleStop)
    ui.clearRunButton.clicked.connect(handleClearRun)
    ui.clearButton.clicked.connect(handleClear)


    t1 = threading.Thread(target=runChecks)
    t1.start()


    window.show()
    ret = app.exec_()
    global converter
    global friend_server
    global run_checks
    run_checks = False
    if converter != None:
        converter.stopExe()
    if friend_server != None:
        friend_server.stopExe()
    t1.join()
    sys.exit(ret)
