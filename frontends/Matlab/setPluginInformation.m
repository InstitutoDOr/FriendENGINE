% sending plugin information
function setPluginInformation(friendObj)

% sending the PLUGIN command and parameters
fprintf(friendObj.mainThread, 'PLUGIN');
fprintf(friendObj.mainThread, 'libROI');
fprintf(friendObj.mainThread, 'no');
fprintf(friendObj.mainThread, 'processROI');
fprintf(friendObj.mainThread, 'initializeROIProcessing');
fprintf(friendObj.mainThread, 'finalizeProcessing');
fprintf(friendObj.mainThread, 'no');
fprintf(friendObj.mainThread, 'no');
% getting the acknowledge
response=fgetl(friendObj.mainThread);
