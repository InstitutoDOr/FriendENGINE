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
 
# changing the MNI mask
engine.setVariable('MNIMask', 'studydirhmat_spm_final.nii');   

# changing the MNI template
engine.setVariable('MNITemplate', 'studydirMNI152_T1_1mm_brain.nii.gz');

# changing the directory of the volumes
engine.setVariable('Prefix', 'outputdirRUN0' + option + os.path.sep + 'DRIN-');   

# changing the activation level
engine.setVariable('ActivationLevel', '0.01');   

# changing the current suffix
engine.setVariable('CurrentRunSuffix', 'RUN0' + option);   
   
# initiating processing
engine.doTrain = False;
engine.doGLM = False;
engine.doFeatureSelection = False;
engine.additionalFeedbacks = 1;

engine.startTheEngine(3, True);
