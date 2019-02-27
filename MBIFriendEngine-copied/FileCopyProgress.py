import sys, os, datetime
from PyQt4 import QtGui, QtCore
import shutil
import os.path


class FileCopyProgress( QtGui.QDialog ):
    '''Custom shutil file copy with progress'''
    def __init__(self, app, parent=None, src=None, dest=None):
        super(FileCopyProgress, self).__init__()
    
        self.src = src
        self.dest = dest
        self.app = app
        self.build_ui()



    def build_ui(self):

        hbox = QtGui.QVBoxLayout()

        lbl_src = QtGui.QLabel('Source: ' + self.src)
        lbl_dest = QtGui.QLabel('Destination: ' + self.dest)
        self.pb = QtGui.QProgressBar()

        self.pb.setMinimum(0)
        self.pb.setMaximum(100)
        self.pb.setValue(0)


        hbox.addWidget(lbl_src)
        hbox.addWidget(lbl_dest)
        hbox.addWidget(self.pb)
        self.setLayout(hbox)
        
        self.setWindowTitle('File copy')

        self.auto_start_timer = QtCore.QTimer.singleShot(2000,lambda : self.copydirs( src=self.src, dst=self.dest, callback_progress=self.progress, callback_copydone=self.copydone ))

        self.exec_()


    def progress(self, copied, tot):
                
        percentage = int((float(copied)/float(tot))*100)
        try:
            self.pb.setValue( percentage )
        except:
            pass

        self.app.processEvents()


    def copydone(self):
        self.pb.setValue( 100 )
        self.close()


    def copydirs(self, src, dst, callback_progress, callback_copydone):
        dirs = []
        files = []
        for dirName, subdirList, fileList in os.walk(src):
            dst_dir = dirName.replace(src, dst+os.path.sep)
            dirs.append(dst_dir)
            for subdir in subdirList:
                if subdir not in dirs:
                    dirs.append(os.path.join(dst_dir,subdir))

            for ff in fileList:
                files.append((os.path.join(dirName,ff),os.path.join(dst_dir,ff)))
        
        ## print len

        print len(files)

        tot = len(files) + len(dirs)

        copied = 0
        for dd in dirs:
            ## replace the source with dest
            
            if not os.access(dd,os.F_OK):
                os.mkdir(dd)
                callback_progress(copied,tot)
            copied = copied +1
        for ff in files:
            shutil.copy(ff[0],ff[1])
            callback_progress(copied,tot)
            copied = copied +1

        callback_copydone()


