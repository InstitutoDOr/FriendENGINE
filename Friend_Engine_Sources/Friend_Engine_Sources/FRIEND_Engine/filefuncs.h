#ifndef FILEFUNCS_H
#define FILEFUNCS_H
#include "defs.h"
#include <vector>
#include <string>
#include <string.h>

using namespace std;

// returns the list of files in dirSource
void listFiles(char *dirSource, vector<string>&entries);

// returns the list of directories in dirSource
void listDirectory(char *dirSource, vector<string>&entries, int fullpath = 1);

// returns the first entry of a directory
string returnFirstFile(char *dirSource);

// returns the first element of a list that has the given substring 
string searchSubstringInList(vector<string>&entries, char *substring, char *noSearchSubString = NULL);

// merge to files in binary mode
void mergeFiles(const char *fromfileA, const char *fromfileB, const char *tofile);

// just copies a file
void copyfile(char *fromfile, char *tofile);

// another version crossplatform
void copyFile(const char *fromfile, const char *tofile);

// just deleting a directory
int removeDirectory(const char *dirName);

// copy a directory
int copyDirectory(const char *dirSource, const char *dirDest);

// verifying if a file is ready to be read
bool isReadable(char *fileName);

// checks if fileName or fileName.gz exists and return the one that exists
int returnFileNameExists(char *fileName, char *fileNameExists);

// checks if a file exists
int fileExists(char *arquivo);

// file Size
int fileSize(char *filename);

bool isReadableSize(char *fileName, int size);

// returns the path of `fileName`. Note `output` must have at least the same size than `fileName`
void extractFilePath(char *file, char *saida);

// window rename directory
#ifdef WINDOWS
int renameDirectory(char *fileName, char *destFileName);
#endif

// move a file respecting the extension
int moveFile(char *fileName, char *destFilename);

// returns the filename part of `fileName`. Note output must have at least the same size than fileName
void extractFileName(char *file, char *saida);

// returns the filename part of `fileName`
char *extractFileName(char *file);

// just changing the extension of the file. Note output must have the necessary size
void changeFileExt(char *Filename, const char* newExt, char *output);

// including a last path delimiter
void includeTrailingPathDelimiter(char *file);

// performs a shell like expansion of fileName
void expandFilename(char *fileName);

#endif
