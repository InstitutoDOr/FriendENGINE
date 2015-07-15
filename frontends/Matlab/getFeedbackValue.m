% get the feedback value. The feedback value was previously calculated and
% stored in a session variable in Friend Engine. This just opens another
% thread, referencing the session to get these values.
function [class, percentage] = getFeedbackValue(friendObj)
    fopen(friendObj.responseThread);
    fprintf(friendObj.responseThread, 'SESSION');
    fprintf(friendObj.responseThread, '%s', friendObj.sessionID);
    response=fgetl(friendObj.responseThread);
    fprintf(friendObj.responseThread, 'TEST'); 
    fprintf(friendObj.responseThread, '%d', friendObj.actualVolume);
    class=str2double(fgetl(friendObj.responseThread));
    percentage=str2double(fgetl(friendObj.responseThread));
    response=fgetl(friendObj.responseThread);
    fclose(friendObj.responseThread);