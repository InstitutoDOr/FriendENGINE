#include <stdlib.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include <stddef.h>
#include <float.h>
//#include <unistd.h>
#include <time.h>  // clock_t, clock, CLOCKS_PER_SEC
#include <stdio.h>
#include "nii_dicom_batch.h"
#include "nii_dicom.h"
#include <math.h>


void defaultOptions(TDCMopts &opts, char *outdir);
int singleDICOM(struct TDCMopts* opts, char *fname);
void transformDirDicom(char *indir, char *output);
void transformDicom(char *filename, char *outdir);
void transformDicom(char *filename);
void transformDicom2(char *filename, char *output, int processDir = 0);
void dropFilenameFromPath(char *path);
int nii_saveCrop(char * niiFilename, struct nifti_1_header hdr, unsigned char* im, struct TDCMopts opts);
struct TDICOMdata readDICOM(char * fname);