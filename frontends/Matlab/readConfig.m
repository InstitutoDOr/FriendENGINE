function readConfig(friendObj, file)
   fileID = fopen(file);                                   % File variable is the filename
   A = fread(fileID, '*char')';                            
   fclose(fileID);                                         % Closes the file
   fopen(friendObj.responseThread);
   fprintf(friendObj.mainThread, 'READCONFIG');
   fprintf(friendObj.responseThread, '%s', size(A, 1)) % Sends length of the file as string 
   fprintf(friendObj.mainThread, '%s', A);              % Sends the content of the file   
   fgetl(friendObj.mainThread);                            % Reading the acknowledge info.
   fclose(friendObj.responseThread);
end