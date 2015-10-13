function friendObj = initializeFriendObject(host, port)
    % graphs variables
    friendObj.rotationx = [];
    friendObj.rotationy = [];
    friendObj.rotationz = [];
    
    friendObj.translationx = [];
    friendObj.translationy = [];
    friendObj.translationz = [];
    friendObj.feedbackValues = [];
    friendObj.rms = [];
    friendObj.percentage = 0;
    friendObj.class = 0.0;
    friendObj.actualVolume = 0;
    friendObj.phase = 0;
    friendObj.pipelineType = 1;
    friendObj.feedbackRun = 1;
    friendObj.additionalFeedbacks = 0;
    
    friendObj.doGLM = 1;
    friendObj.doFEATURESELECTION = 1;
    friendObj.doTRAIN = 0;
    
    % comunication protocol variables
    friendObj.mainThread = 0;
    friendObj.responseThread = 0;
    friendObj.sessionID = '';
    
    timeOut = 1;
    friendObj.mainThread = tcpip(host, port);
    set(friendObj.mainThread, 'TimeOut', timeOut);
	friendObj.mainThread.OutputBufferSize = 10000; % may need to be larger

    friendObj.responseThread = tcpip(host, port);
    set(friendObj.responseThread, 'TimeOut', timeOut);    
end
