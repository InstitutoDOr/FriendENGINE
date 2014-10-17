import socket;
import sys;
import time;
import datetime;

#communication variables
HOST = '127.0.0.1';
PORT = 5678;
TR = 2.0;
sessionID = '';
PIPELINE = 1;
Timeout = 0.3; # I have used 0.05 in the same operational system (Windows). Tested in Mac and Linux (Cent OS), running the engine in Windows.
bufferLines=tuple();

# time stamp
def timestamp():
   ts = time.time();
   st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S');
   return st;

#function to read from a socket
def readsocket(socket):
   global bufferLines;
   if (len(bufferLines)==0): 
      totalData='';
      data='';
      try:
         socket.settimeout(Timeout);
         data = socket.recv(1024);
         while (data):
            totalData +=data;
            data = socket.recv(1024);
      except:
         pass;
      bufferLines=totalData.split('\n');
      if (len(bufferLines) > 0):
         bufferLines.pop(-1);
   if (len(bufferLines) > 0):
      result=bufferLines.pop(0);
   else:
      result = '';
   return result;
   
# function to create a session
def createSession(mainThread):
   mainThread.send('NEWSESSION\n')
   sessId = readsocket(mainThread);
   response = readsocket(mainThread);
   return sessId;

# function to set the plugin information
def setPlugInInformation(mainThread):
# sending the PLUGIN command and parameters
   mainThread.send('PLUGIN\n');
   mainThread.send('libBrainDecoding\n');
   mainThread.send('trainSVM\n');
   mainThread.send('testSVM\n');
   mainThread.send('initSVM\n');
   mainThread.send('finalSVM\n');
   mainThread.send('no\n');
   mainThread.send('no\n');

# getting the acknowledge
   response=readsocket(mainThread);

# get the status of a command   
def getResponse(*arg):
   responseThread = socket.socket(socket.AF_INET, socket.SOCK_STREAM);
   responseThread.connect((HOST, PORT));
   responseThread.send('SESSION\n');
   responseThread.send(sessionID + '\n');

   response=readsocket(responseThread);
   if (response == 'OK'):
      responseThread.send(arg[0]);
      if (len(arg) > 1):
         responseThread.send(str(arg[1]) + '\n'); 
   try:
      response=readsocket(responseThread);
      if (len(arg) > 1):
	     # this is the true acknowledge 
         ack=readsocket(responseThread);
   except:
      print("Time out !!!");
   responseThread.close();
   return response;
   
# get feedback information
def getFeedbackValue(actualVolume):
    global bufferLines;
    responseThread = socket.socket(socket.AF_INET, socket.SOCK_STREAM);
    responseThread.connect((HOST, PORT));
    print("Sending TEST command");
    responseThread.send('SESSION\n');
    responseThread.send(sessionID + '\n');

    response=readsocket(responseThread);
    responseThread.send('TEST\n');
    responseThread.send(str(actualVolume) + '\n');

    classe=readsocket(responseThread);
    percentage=readsocket(responseThread);
	# reading acknowledge
    response=readsocket(responseThread);
    return (classe, percentage);
	
def processPhase(phase, actualVolume, mainThread):
   global bufferLines;
   if (phase == 1):
      # sending PREPROC command
      print("sending PREPROC command");
      mainThread.send('NBPREPROC\n');
      response=readsocket(mainThread);
      phase = 15;

   if (phase == 15):
      # sees if PREPROC terminates
      print("Querying PREPROC termination status.");
      response=getResponse('PREPROC\n');
      print("Response received = %s" % (response));
      if (response == 'OK'):
         phase = 2;
 
   if (phase == 2):
      # sending FEEDBACK command
      print("sending FEEDBACK command");
      if (PIPELINE==1):	  
          mainThread.send('NBPIPELINE\n');
      else:
          mainThread.send('NBFEEDBACK\n');
      response=readsocket(mainThread);
      print("Response received = %s" % (response));
      actualVolume=1;
      phase = 25;

   if (phase == 25):
      #retrieving Graph Parameters
      print("Querying graph parameters of volume %d." % (actualVolume));
      response=getResponse('GRAPHPARS\n', actualVolume);
      values=response.split(';');
      if (response == 'END'):
         if (PIPELINE==1):
            # GLM session. This is the blocked version. Its sends the response at the end of the processing. 
            # NBGLM is the non-blocked version
            mainThread.setblocking(1);
            mainThread.send('GLM\n');
            response=mainThread.recv(3);

            # feature selection. Non-blocked version is NBFEATURESELECTION
            mainThread.send('FEATURESELECTION\n');
            response=mainThread.recv(3);
	 
            # calling the plug-in train function, which builds the svm model
            mainThread.send('TRAIN\n');
            response=mainThread.recv(3);

         # ending the session
         mainThread.send('ENDSESSION\n');
         mainThread.send(sessionID+'\n');
         response=readsocket(mainThread);
         phase = 100;
      elif (len(values) == 9):
         print("%s : Volume = %s"  % (timestamp(), values[1]));
         print("ROT   (X = %s,  Y = %s, Z = %s)"  % (values[2], values[3], values[4]));
         print("TRANS (X = %s,  Y = %s, Z = %s)"  % (values[5], values[6], values[7]));
         print("RMS      = %s"  % (values[8]));
         if (PIPELINE == 2):
            classe, feedback = getFeedbackValue(actualVolume);
            print("Class returned %s  Feedback=%s\n" % (classe, feedback));
         actualVolume = actualVolume + 1;
      else:
         print("Response received = %s" % (response));
   return (phase, actualVolume); 


#socket variables
mainThread = socket.socket(socket.AF_INET, socket.SOCK_STREAM);

#connection main socket
mainThread.connect((HOST, PORT));

# creating a session and getting sessionID
sessionID = createSession(mainThread);
if (PIPELINE==2):
   # changing the mask type
   mainThread.send('SET\n');
   mainThread.send('ModelRunSuffix\n');
   mainThread.send('RUN01\n');
   response=readsocket(mainThread);
   
   # changing the mask
   mainThread.send('SET\n');
   mainThread.send('CurrentRunSuffix\n');
   mainThread.send('RUN02\n');
   response=readsocket(mainThread);
   
   # changing the directory of the volumes
   mainThread.send('SET\n');
   mainThread.send('Prefix\n');
   mainThread.send('outputdirRUN02/DRIN-\n');
   response=readsocket(mainThread);
   
print("sessionID received = %s" % (sessionID));

# configuring plugIn information
setPlugInInformation(mainThread);

# initiating processing
phase = 1;
actualVolume = 1;
while (phase < 100):
   phase, actualVolume = processPhase(phase, actualVolume, mainThread);
   time.sleep(TR / 7.0);
