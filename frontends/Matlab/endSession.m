function friendObj = endSession(friendObj)
   fprintf(friendObj.mainThread, 'ENDSESSION');
   fprintf(friendObj.mainThread, '%s', friendObj.sessionID);
   fgetl(friendObj.mainThread);
   fclose(friendObj.mainThread);
   friendObj.phase = 100;
end;