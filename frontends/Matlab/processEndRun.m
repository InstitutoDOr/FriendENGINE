function friendObj = processEndRun(friendObj)
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