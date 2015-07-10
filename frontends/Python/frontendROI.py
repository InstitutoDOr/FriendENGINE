import socket;
import sys;
import os;
from friendEngine import Engine;

engine = Engine();

# creating the main socket variable
engine.connectEngine();

# creating a session and getting sessionID
engine.createSession();
print("sessionID received = %s" % (engine.sessionID));

# reading the actual run number
if (len(sys.argv)>1):
   option = str(sys.argv[1]);
else:
   option = "1";

PIPELINE = int(option);
 
if (PIPELINE > 1):
   # changing the mask type
   engine.setVariable('ActivationLevelMaskType', '1');   
   # changing the mask
   engine.setVariable('ActivationLevelMask', 'glmdirtstats_features_RUN01_bin');   
   # changing the directory of the volumes
   engine.setVariable('Prefix', 'outputdirRUN0' + option + os.path.sep + 'DRIN-');   
   # changing the current suffix
   engine.setVariable('CurrentRunSuffix', 'RUN0' + option);   
   
# initiating processing
engine.startTheEngine(1, PIPELINE!=1);