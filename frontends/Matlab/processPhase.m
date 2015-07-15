% Executes a phase in the Friend Engine Pipeline : (PREPROC, FEEDBACK,
% query for motion correction parameters and results
function friendObj = processPhase(friendObj)
if (friendObj.phase == 1)
    % sending PREPROC command
    fprintf(friendObj.mainThread, 'NBPREPROC');
    fgetl(friendObj.mainThread);
    friendObj.phase = 15;
end

if (friendObj.phase == 15)
    % sees if PREPROC terminates
    response=getResponse(friendObj, 'PREPROC');
    if (strcmp(response, 'OK') == 1)
       friendObj.phase = 2;
    end;
end

if (friendObj.phase == 2)
    % sending FEEDBACK command
    if (friendObj.feedbackRun == 1)
       fprintf(friendObj.mainThread, 'NBFEEDBACK');
    else
       fprintf(friendObj.mainThread, 'NBPIPELINE');
    end;
    response=fgetl(friendObj.mainThread);
    friendObj.actualVolume=1;
    friendObj.phase = 25;
end

if (friendObj.phase == 25)
    % retrieving Graph Parameters
    response=getResponse(friendObj, 'GRAPHPARS', friendObj.actualVolume);
    if (strcmp(response, 'END') == 1)
       friendObj=processEndRun(friendObj);
       return;
    end;
    
    if (strcmp(response, 'NOK.') == 0)
       tokens=regexp(response, ';', 'split');
       if (size(tokens, 2) == 9)

            % getting the feedback
            [friendObj.class, friendObj.percentage] = getFeedbackValue(friendObj);

            % incrementing the to be processed volume
            friendObj.actualVolume = friendObj.actualVolume + 1;

            % updating graph vars
            friendObj.rotationx = [friendObj.rotationx str2double(tokens{3})];
            friendObj.rotationy = [friendObj.rotationy str2double(tokens{4})];
            friendObj.rotationz = [friendObj.rotationz str2double(tokens{5})];

            friendObj.translationx = [friendObj.translationx str2double(tokens{6})];
            friendObj.translationy = [friendObj.translationy str2double(tokens{7})];
            friendObj.translationz = [friendObj.translationz str2double(tokens{8})];

            friendObj.rms = [friendObj.rms str2double(tokens{9})];
        end;
    end;
end
