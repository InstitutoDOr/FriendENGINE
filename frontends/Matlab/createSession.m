% create session
function sessionID = createSession(friendObj)
fprintf(friendObj.mainThread, 'NEWSESSION');
% reading session id
sessionID = fgetl(friendObj.mainThread);
friendObj.sessionID = sessionID;
% reading ok
response=fgetl(friendObj.mainThread);