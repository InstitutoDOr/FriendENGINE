function friendObj = stopSession(friendObj)
   fopen(friendObj.responseThread);
   fprintf(friendObj.responseThread, 'STOPSESSION');
   fprintf(friendObj.responseThread, '%s', friendObj.sessionID);
   fgetl(friendObj.responseThread);
   fclose(friendObj.responseThread);
   friendObj.doGLM = 0;
   friendObj.doFEATURESELECTION = 0;
   friendObj.doTRAIN = 0;   
end