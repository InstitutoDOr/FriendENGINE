function configureVariable(friendObj, variable, value)
   fprintf(friendObj.mainThread, 'SET');
   fprintf(friendObj.mainThread, variable);
   fprintf(friendObj.mainThread, value);
   fgetl(friendObj.mainThread);
end