import threading
from PyQt4 import QtCore
from PyQt4.QtCore import QObject
import paramiko
import os
import os.path
import time
import string

class SignalHandler(QObject):

    updateInfoSig = QtCore.pyqtSignal('QString',bool)

    def __init__(self):
        super(SignalHandler,self).__init__()

    def updateInfo(self,data,err=False):
        self.updateInfoSig.emit(data,err)
    

class RemoteRunner(threading.Thread):

   
    def __init__(self,cmd, args, keyfName):
        super(RemoteRunner,self).__init__()
        self.cmd = cmd
        self.args = args
        self.key = keyfName
        self.sig_handler = SignalHandler()
        self.running = False


    def stopExe(self):
        self.running =False

    def run(self):

        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        #client.load_host_keys(os.path.expanduser('~/.ssh/known_hosts'))
        k = paramiko.RSAKey.from_private_key_file(self.key)
        client.connect("carbon2.mbi.monash.edu", username="nfb", pkey=k)
        stdin, stdout, stderr = client.exec_command(self.cmd+" "+" ".join(self.args),get_pty=True)
        print "Running"
        self.running = True
        chan = stdout.channel

        while self.running and not chan.exit_status_ready():
            if stdout.channel.recv_ready():
                print "Recv Out"
                data = stdout.channel.recv(1024)
                self.sig_handler.updateInfo(data,False)
                #    data = stdout.channel.recv(1024)
            if stderr.channel.recv_ready():
                print "Recv Err"
                data = stderr.channel.recv(1024)
                self.sig_handler.updateInfo(data,True)
                 #   data = stdout.channel.recv(1024)
         #   print self.running
            if not self.running:
                break
            time.sleep(0.01)

        print "Exit"
        client.close()
        print "Closing client " + self.cmd
