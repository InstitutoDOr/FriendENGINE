% get a process terminated ? response
function [response] = getResponse(friendObj, varargin)

fopen(friendObj.responseThread);
fprintf(friendObj.responseThread, 'SESSION');
fprintf(friendObj.responseThread, '%s', friendObj.sessionID);
response=fgetl(friendObj.responseThread);
if (strcmp(response, 'OK') == 1)
   fprintf(friendObj.responseThread, '%s', varargin{1});
   if (length(varargin) > 1)
      fprintf(friendObj.responseThread, '%d', varargin{2}); 
   end;
   response = fgetl(friendObj.responseThread);
   if (length(varargin) > 1)
      fgetl(friendObj.responseThread);      
   end;
end;
fclose(friendObj.responseThread);