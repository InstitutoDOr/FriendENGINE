#include "zlib.h"

#ifndef GZ_SUFFIX
#  define GZ_SUFFIX ".gz"
#endif
#define SUFFIX_LEN (sizeof(GZ_SUFFIX)-1)

#define BUFLEN      16384
#define MAX_NAME_LEN 1024

#ifdef MAXSEG_64K
#  define local static
/* Needed for systems with limitation on stack size. */
#else
#  define local
#endif


void gz_compress(FILE *in, gzFile out);
void gz_uncompress(gzFile in, FILE *out);
void file_compress(char  *file, char  *mode);
int file_uncompress(char  *file, char *outputFile);