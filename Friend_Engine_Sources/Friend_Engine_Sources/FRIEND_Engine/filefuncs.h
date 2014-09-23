#ifndef FILEFUNCS_H
#define FILEFUNCS_H
#include "defs.h"

// just copies a file
void copyfile(char *fromfile, char *tofile);

// just deleting a directory
int removeDirectory(const char *dirName);

// verifying if a file is ready to be read
bool isReadable(char *fileName);

// checks if fileName or fileName.gz exists and return the one that exists
int returnFileNameExists(char *fileName, char *fileNameExists);

// checks if a file exists
int fileExists(char *arquivo);

// returns the path of `fileName`. Note `output` must have at least the same size than `fileName`
void extractFilePath(char *file, char *saida);

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
