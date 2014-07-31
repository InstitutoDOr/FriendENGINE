#ifndef WINDOWS
#include <dirent.h>
#ifdef UNIX
#ifndef LINUX
#include "copyfile.h"
#endif
#endif
#include <wordexp.h>
#else
#include "dirent.h"
#include "io.h"
#endif

#include "filefuncs.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

using namespace std;

// performs a shell like expansion of fileName
void expandFilename(char *fileName)
{
#ifndef WINDOWS
   wordexp_t exp_result;
   wordexp(fileName, &exp_result, 0);
   fileName[0]=0;
   strcpy(fileName, exp_result.we_wordv[0]);
   wordfree(&exp_result);
#endif
}

// just copying the files
void copyfile(char *fromfile, char *tofile)
{
#ifdef UNIX
#ifndef LINUX
   copyfile_state_t state;
   state = copyfile_state_alloc();
   copyfile(fromfile, tofile, state, COPYFILE_DATA | COPYFILE_XATTR);
   copyfile_state_free(state);
#endif
#endif
}

// just deleting a directory
int removeDirectory(const char *dirName)
{
   DIR *dir;
   struct dirent *entry;
   char path[PATH_MAX];
   
   dir = opendir(dirName);
   if (dir == NULL) {
      return 0;
   }
   
   while ((entry = readdir(dir)) != NULL) {
      if (strcmp(entry->d_name, ".") && strcmp(entry->d_name, "..")) {
         sprintf(path, "%s/%s", dirName, entry->d_name);
         if (entry->d_type == DT_DIR) {
            removeDirectory(path);
         }        
         remove(path);
      }
      
   }
   closedir(dir);
   remove(dirName);
   
   return 1;
}

// verifying if a file is ready to be read
bool isReadable(char *fileName)
{
   fstream filen;
   int r;
   filen.open(fileName, ios::in | ios::out);
   if (filen.good()) r =true;
   else r = false;
   filen.close();
   return r;
}

// checks if fileName or fileName.gz exists and return the one that exists
int returnFileNameExists(char *fileName, char *fileNameExists)
{
   int resp;
#ifdef WIN32
   if (_access(fileName, 0) != -1)
#else
      if (access(fileName, 0) != -1)
#endif
         resp=1;
      else resp=0;
   
   if (resp==0)
   {
      char auxfilename[2048];
      sprintf(auxfilename, "%s%s", fileName, ".gz");
#ifdef WIN32
      if (_access(auxfilename, 0) != -1)
#else
         if (access(auxfilename, 0) != -1)
#endif
            resp=1;
         else resp=0;
      if (resp)
         sprintf(fileNameExists, "%s", auxfilename);
   }
   else sprintf(fileNameExists, "%s", fileName);
   return resp;
}

// checks if a file exists
int fileExists(char *fileName)
{
   int resp;
   if (fileName==NULL) return 0;
	#ifdef WIN32
   if (_access(fileName, 0) != -1)
   #else
   if (access(fileName, 0) != -1)
   #endif
      resp=1;
	else resp=0;
   
   if (resp==0)
   {
      char gzFileName[2048];
      sprintf(gzFileName, "%s%s", fileName, ".gz");
      #ifdef WIN32
      if (_access(gzFileName, 0) != -1)
      #else
      if (access(gzFileName, 0) != -1)
      #endif
         resp=1;
      else resp=0;
      if (resp)
         strcpy(fileName, gzFileName);
   }
   return resp;
}

// returns the path of `fileName` without the last PATHSEPARATOR. Note `output` must have at least the same size than `fileName`
void extractFilePath(char *fileName, char *output)
{
	strcpy(output, fileName);
	int i;
	for(i=strlen(output)-1;i >= 0; i--)
		if (output[i] == PATHSEPCHAR)
		{
			output[i]=0;
			break;
		}
	if (i <= 0) output[0] = 0;
}

// returns the filename part of `fileName`. Note `output` must have at least the same size than `fileName`
void extractFileName(char *fileName, char *output)
{
	for(int i=strlen(fileName)-1; i >= 0; i--)
		if (fileName[i] == PATHSEPCHAR)
		{
			strcpy(output, (fileName + i+1));
			break;
		}
}

// returns the path of `fileName`. Please free the memory. I cant garantee if I do
char *extractFileName(char *file)
{
   char *fileName = (char *) malloc((strlen(file)+1) * sizeof(char));
   for (int i=strlen(file)-1;i >= 0; i--)
       if (file[i] == PATHSEPCHAR) 
       {
          strcpy(fileName, (file + i+1));
          break;
       }
   return fileName;
}

// just changing the extension of the file. Note output must have the necessary size
void changeFileExt(char *filename, const char* newExt, char *output)
{
    strcpy(output, filename);
    char *end=strrchr(output, '.');
    if (end)
    {
        end[0]=0;
        strcat(output, newExt);
    }
}

// including a last path delimiter
void includeTrailingPathDelimiter(char *path)
{
   if (path[strlen(path)-1] != PATHSEPCHAR)
   {
      int pathSize=strlen(path);
      path[pathSize] = PATHSEPCHAR;
      path[pathSize+1] = 0;
   }
}
