import socket;
import time;
import datetime;

class Engine(object):
   def __init__(self):

      #communication variables
      self.Timeout = 0.1; 

      # processing variables   
      self.TR = 2.0;
      self.sessionID = '';
      self.actualVolume=1;
      self.feedbackRun=0;
      self.phase=1;
      self.bufferLines=tuple();
      self.doTrain = True;
      self.doGLM = True;
      self.doFeatureSelection = True;
      self.additionalFeedbacks = 0;
   
   # time stamp
   def timestamp(self):
      ts = time.time();
      st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S');
      return st;

   #function to read from a socket
   def readsocket(self, socket):
      if (len(self.bufferLines)==0): 
         totalData='';
         data='';
         try:
            socket.settimeout(self.Timeout);
            data = socket.recv(1024);
            while (data):
               totalData +=data;
               data = socket.recv(1024);
         except:
            pass;
         self.bufferLines=totalData.split('\n');
         if (len(self.bufferLines) > 0):
            self.bufferLines.pop(-1);
      if (len(self.bufferLines) > 0):
         result=self.bufferLines.pop(0);
      else:
         result = '';
      return result;

   def connectEngine(self, engineIp = '127.0.0.1', portNumber = 5678):
      #Creating the socket variable
      self.mainThread = socket.socket(socket.AF_INET, socket.SOCK_STREAM);

      #connection main socket
      self.connectionInfo = (engineIp, portNumber);
      self.mainThread.connect(self.connectionInfo);
      
   # function to create a session
   def createSession(self):
      self.mainThread.send('NEWSESSION\n')
      self.sessionID = self.readsocket(self.mainThread);
      response = self.readsocket(self.mainThread);
      return self.sessionID;

   # function to set the plug-in information
   def setROIPlugIn(self):
      # sending the PLUGIN command and parameters
      self.mainThread.send('PLUGIN\n');
      self.mainThread.send('libROI\n');
      self.mainThread.send('no\n');
      self.mainThread.send('processROI\n');
      self.mainThread.send('initializeROIProcessing\n');
      self.mainThread.send('finalizeProcessing\n');
      self.mainThread.send('no\n');
      self.mainThread.send('no\n');

      # getting the acknowledge
      response = self.readsocket(self.mainThread);
      return response;   
   
   # function to set the plug-in information
   def setConnectivityPlugIn(self):
      # sending the PLUGIN command and parameters
      self.mainThread.send('PLUGIN\n');
      self.mainThread.send('libConnectivity\n');
      self.mainThread.send('buildROIs\n');
      self.mainThread.send('calculateFeedback\n');
      self.mainThread.send('initializeFunctionalConectivity\n');
      self.mainThread.send('finalizeFunctionalConectivity\n');
      self.mainThread.send('no\n');
      self.mainThread.send('no\n');

      # getting the acknowledge
      response = self.readsocket(self.mainThread);
      return response;   
	  
   # function to set the plug-in information
   def setSVMPlugIn(self):
      # sending the PLUGIN command and parameters
      self.mainThread.send('PLUGIN\n');
      self.mainThread.send('libBrainDecoding\n');
      self.mainThread.send('trainSVM\n');
      self.mainThread.send('testSVM\n');
      self.mainThread.send('initSVM\n');
      self.mainThread.send('finalSVM\n');
      self.mainThread.send('no\n');
      self.mainThread.send('no\n');

      # getting the acknowledge
      response = self.readsocket(self.mainThread);
      return response;   
   
   # function to set the plug-in information
   def setMotorPlugIn(self):
      # sending the PLUGIN command and parameters
      self.mainThread.send('PLUGIN\n');
      self.mainThread.send('libMotor\n');
      self.mainThread.send('no\n');
      self.mainThread.send('processMotorROI\n');
      self.mainThread.send('initializeMotorProcessing\n');
      self.mainThread.send('finalizeMotorProcessing\n');
      self.mainThread.send('no\n');
      self.mainThread.send('no\n');

      # getting the acknowledge
      response = self.readsocket(self.mainThread);
      return response;   


   # function to set the plug-in information
   def setPlugInInformation(self, plugInType):
      if (plugInType==1):
         return self.setROIPlugIn();
	  
      if (plugInType==2):
         return self.setSVMPlugIn();
	  
      if (plugInType==3):
         return self.setMotorPlugIn();

      if (plugInType==4):
         return self.setConnectivityPlugIn();

   # get the status of a command   
   def getResponse(self, *arg):
      responseThread = socket.socket(socket.AF_INET, socket.SOCK_STREAM);
      responseThread.connect(self.connectionInfo);
      responseThread.send('SESSION\n');
      responseThread.send(self.sessionID + '\n');

      response = self.readsocket(responseThread);
      if (response == 'OK'):
         responseThread.send(arg[0]);
         if (len(arg) > 1):
            responseThread.send(str(arg[1]) + '\n'); 
      try:
         response = self.readsocket(responseThread);
         if (len(arg) > 1):
	        # this is the true acknowledge 
            ack = self.readsocket(responseThread);
      except:
         print("Time out !!!");
      responseThread.close();
      return response;
   
   # get feedback information
   def getFeedbackValue(self):
      global bufferLines;
      responseThread = socket.socket(socket.AF_INET, socket.SOCK_STREAM);
      responseThread.connect(self.connectionInfo);
      print("Sending TEST command");
      responseThread.send('SESSION\n');
      responseThread.send(self.sessionID + '\n');

      response = self.readsocket(responseThread);
      responseThread.send('TEST\n');
      responseThread.send(str(self.actualVolume) + '\n');

      classe = self.readsocket(responseThread);
      percentage = self.readsocket(responseThread);
      if (self.additionalFeedbacks > 0):
         for x in xrange(0, self.additionalFeedbacks):
             response = self.readsocket(responseThread);             
	  # reading acknowledge
      response = self.readsocket(responseThread);
      responseThread.close();
      return (classe, percentage);

   def sendCommand(self, command, timeout=3):
      self.mainThread.send(command + '\n');
      response = self.mainThread.recv(timeout);
      return response;
   
   def endSession(self):
      self.mainThread.send('ENDSESSION\n');
      self.mainThread.send(self.sessionID+'\n');
      response = self.readsocket(self.mainThread);
      return response;
   
   def stopSession(self):
      responseThread = socket.socket(socket.AF_INET, socket.SOCK_STREAM);
      responseThread.connect(self.connectionInfo);
      responseThread.send('STOPSESSION\n');
      responseThread.send(self.sessionID+'\n');
      response = self.readsocket(responseThread);
      responseThread.close();
      self.doTrain = False;
      self.doGLM = False;
      self.doFeatureSelection = False;
      return response;
	  
   def printMotionParameters(self, motionparams):
      if (len(motionparams) == 9):
         print("%s : Volume = %s"  % (self.timestamp(), motionparams[1]));
         print("ROT   (X = %s,  Y = %s, Z = %s)"  % (motionparams[2], motionparams[3], motionparams[4]));
         print("TRANS (X = %s,  Y = %s, Z = %s)"  % (motionparams[5], motionparams[6], motionparams[7]));
         print("RMS      = %s"  % (motionparams[8]));


   def setVariable(self, variable, value):
      self.mainThread.send('SET\n');
      self.mainThread.send(variable + '\n');
      self.mainThread.send(value + '\n');
      response = self.readsocket(self.mainThread);
      return response;

   def processEndRun(self):
      self.mainThread.setblocking(1);

      if (self.doGLM):
         # GLM session. This is the blocked version. Its sends the response at the end of the processing. 
         # NBGLM is the non-blocked version
         print("Sending GLM command.");
         self.sendCommand('GLM');

      if (self.doFeatureSelection):
         # feature selection. 
         print("Sending feature selection command.");
         self.sendCommand('FEATURESELECTION');
 
      if (self.doTrain):
         # train function. If the plug-in does not implement this, nothing is done.
         print("Sending train command.");
         self.sendCommand('TRAIN');
		
      # ending the session
      self.endSession();

   
   def processPhase(self, feedbackRun, oneStep=False):
      if (self.phase == 1):
         # sending PREPROC command
         # print("sending PREPROC command");
         self.mainThread.send('NBPREPROC\n');
         response = self.readsocket(self.mainThread);
         self.phase = 15;
         if (oneStep):
            return (self.phase, self.actualVolume)

      if (self.phase == 15):
         # sees if PREPROC terminates
         # print("Querying PREPROC termination status.");
         response=self.getResponse('PREPROC\n');
         print("Response received = %s" % (response))
         if (response == 'OK'):
            self.phase = 2;
            if (oneStep):
               return (self.phase, self.actualVolume)
 
      if (self.phase == 2):
         # processing volumes
         print("sending PIPELINE command");
         if (feedbackRun):	  
            self.mainThread.send('NBFEEDBACK\n')
         else:
            self.mainThread.send('NBPIPELINE\n')
         response = self.readsocket(self.mainThread)
         print("Response received = %s" % (response))
         self.actualVolume=1;
         self.phase = 25;
         if (oneStep):
            return (self.phase, self.actualVolume)

      if (self.phase == 25):
         #retrieving Graph Parameters
         print("Querying graph parameters of volume %d." % (self.actualVolume))
         response=self.getResponse('GRAPHPARS\n', self.actualVolume)
         values=response.split(';')
         if (response == 'END'):
            print("END signal received.")
            self.processEndRun();		 
            self.phase = 100;
         elif (len(values) == 9):
            self.printMotionParameters(values);
            if (feedbackRun):
               classe, feedback = self.getFeedbackValue()
               print("Feedback=%s\n" % (feedback))
            self.actualVolume = self.actualVolume + 1
         else:
            print("Response received = %s" % (response))
      return (self.phase, self.actualVolume)

   def startTheEngine(self, plugInType, feedbackRun):
      # configuring plugIn information
      self.setPlugInInformation(plugInType);
   
      self.phase = 1;
      self.actualVolume = 1;
      while (self.phase < 100):
         self.processPhase(feedbackRun);
         time.sleep(self.TR / 7);   