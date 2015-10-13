function friendObj = processEndRun(friendObj)
    % setting timeout to 5 minutes due to complicated processes
    set(friendObj.mainThread, 'TimeOut', 300);
    if (friendObj.doGLM == 1)
       fprintf(friendObj.mainThread, 'GLM');
       fgetl(friendObj.mainThread);
    end;
    
    if (friendObj.doFEATURESELECTION == 1)
      fprintf(friendObj.mainThread, 'FEATURESELECTION');
      fgetl(friendObj.mainThread);
    end;

   if (friendObj.doTRAIN == 1)
      fprintf(friendObj.mainThread, 'TRAIN');
      fgetl(friendObj.mainThread);
   end;
   endSession(friendObj);
end