% get the feedback value. The feedback value was previously calculated and
% stored in a session variable in Friend Engine. This just opens another
% thread, referencing the session to get these values.
function friendObj = getFeedbackValue(friendObj)
    fopen(friendObj.responseThread);
    fprintf(friendObj.responseThread, 'SESSION');
    fprintf(friendObj.responseThread, '%s', friendObj.sessionID);
    response=fgetl(friendObj.responseThread);
    fprintf(friendObj.responseThread, 'TEST');
    fprintf(friendObj.responseThread, '%d', friendObj.actualVolume);
    friendObj.class=str2double(fgetl(friendObj.responseThread));
    
    friendObj.percentage=str2double(fgetl(friendObj.responseThread));
    friendObj.feedbackValues(1) = friendObj.percentage;
    
    for i=1:friendObj.additionalFeedbacks
        friendObj.feedbackValues(i+1) = str2double(fgetl(friendObj.responseThread)); 
    end;
    
    fprintf('Volume = %d\n', friendObj.actualVolume);
    for i=1:(friendObj.additionalFeedbacks+1)
        fprintf('Feedback Value %d = %f\n', i, friendObj.feedbackValues(i)); 
    end;
    fprintf('\n');
    
    response=fgetl(friendObj.responseThread);
    fclose(friendObj.responseThread);
