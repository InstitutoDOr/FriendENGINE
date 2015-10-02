% sending plugin information
function friendObj = setPluginInformation(friendObj, plugInType)
 
if (plugInType == 1)
    % sending the PLUGIN command and parameters
    fprintf(friendObj.mainThread, 'PLUGIN');
    fprintf(friendObj.mainThread, 'libROI');
    fprintf(friendObj.mainThread, 'no');
    fprintf(friendObj.mainThread, 'processROI');
    fprintf(friendObj.mainThread, 'initializeROIProcessing');
    fprintf(friendObj.mainThread, 'finalizeProcessing');
    fprintf(friendObj.mainThread, 'no');
    fprintf(friendObj.mainThread, 'no');
end;

if (plugInType == 2)
    % sending the PLUGIN command and parameters
    fprintf(friendObj.mainThread, 'PLUGIN');
    fprintf(friendObj.mainThread, 'libBrainDecoding');
    fprintf(friendObj.mainThread, 'trainSVM');
    fprintf(friendObj.mainThread, 'testSVM');
    fprintf(friendObj.mainThread, 'initSVM');
    fprintf(friendObj.mainThread, 'finalSVM');
    fprintf(friendObj.mainThread, 'no');
    fprintf(friendObj.mainThread, 'no');
end;

if (plugInType == 3)
    % sending the PLUGIN command and parameters
    fprintf(friendObj.mainThread, 'PLUGIN');
    fprintf(friendObj.mainThread, 'libMotor');
    fprintf(friendObj.mainThread, 'no');
    fprintf(friendObj.mainThread, 'processMotorROI');
    fprintf(friendObj.mainThread, 'initializeMotorProcessing');
    fprintf(friendObj.mainThread, 'finalizeMotorProcessing');
    fprintf(friendObj.mainThread, 'no');
    fprintf(friendObj.mainThread, 'no');
    friendObj.additionalFeedbacks = 1;
end;

if (plugInType == 4)
    % sending the PLUGIN command and parameters
    fprintf(friendObj.mainThread, 'PLUGIN');
    fprintf(friendObj.mainThread, 'libConnectivity');
    fprintf(friendObj.mainThread, 'buildROIs');
    fprintf(friendObj.mainThread, 'calculateFeedback');
    fprintf(friendObj.mainThread, 'initializeFunctionalConectivity');
    fprintf(friendObj.mainThread, 'finalizeFunctionalConectivity');
    fprintf(friendObj.mainThread, 'no');
    fprintf(friendObj.mainThread, 'no');
end;
% getting the acknowledge
response=fgetl(friendObj.mainThread);