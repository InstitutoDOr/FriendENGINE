#include "dcm2niiInterface.h"
#include "filefuncs.h"

void defaultOptions(TDCMopts &opts, char *outdir)
{
	strcpy(opts.indir, "");
	strcpy(opts.outdir, "");
	strcpy(opts.imageComments, "");
	opts.isOnlySingleFile = false; //convert all files in a directory, not just a single file
	opts.isRenameNotConvert = false;
	opts.isForceStackSameSeries = false;
	opts.isIgnoreDerivedAnd2D = false;
	opts.isPhilipsFloatNotDisplayScaling = true;
	opts.isCrop = false;
	opts.isGz = false;
	opts.isSave3D = false;
	opts.dirSearchDepth = 5;
	opts.gzLevel = -1;

	opts.isMaximize16BitRange = false; //e.g. if INT16 image has range 0..500 scale to be 0..50000 with hdr.scl_slope =  hdr.scl_slope * 0.01
	opts.isFlipY = true; //false: images in raw DICOM orientation, true: image rows flipped to cartesian coordinates
	opts.isRGBplanar = false; //false for NIfTI (RGBRGB...), true for Analyze (RRR..RGGG..GBBB..B)
	opts.isCreateBIDS = true;
	opts.isOnlyBIDS = false;
	opts.isSortDTIbyBVal = false;
	opts.isAnonymizeBIDS = true;
	opts.isCreateText = false;
	opts.isVerbose = false;
	opts.isTiltCorrect = true;
	opts.numSeries = 0;
	memset(opts.seriesNumber, 0, sizeof(opts.seriesNumber));
	strcpy(opts.filename, "%f_%p_%t_%s");

	opts.isCreateBIDS = false;
	opts.isForceStackSameSeries = false;
	opts.isOnlySingleFile = false;
	opts.isCreateText = false;
	opts.isVerbose = 0;
	opts.isCrop = false;
	opts.pigzname[0] = 0;
	opts.filename[0] = 0;
	opts.isGz = false;
	opts.isIgnoreDerivedAnd2D = false;
	strcpy(opts.outdir, outdir);
}


void transformDicom(char *filename)
{
	char outdir[1024];
	strcpy(outdir, filename);
	dropFilenameFromPath(outdir);
	transformDicom(filename, outdir);
}

void transformDicom(char *filename, char *outdir)
{
	TDCMopts opts;
	char auxFile[1024];

	defaultOptions(opts, outdir);

	changeFileExt(filename, "", auxFile);
	extractFileName(auxFile, opts.filename);
//	opts.isVerbose = 2;
	if (0)
	{
		strcpy(opts.indir, filename);
		nii_loadDir(&opts);
	}
	else singleDICOM(&opts, filename);
}

void transformDirDicom(char *indir, char *output)
{
	TDCMopts opts;
	char auxFile[1024], outdir[1024];

	extractFilePath(output, outdir);
	defaultOptions(opts, outdir);

	changeFileExt(output, "", auxFile);
	extractFileName(auxFile, opts.filename);

	strcpy(opts.indir, indir);
	opts.isGz = 1;
	opts.isCrop = 0;
	opts.isRenameNotConvert = false;
	opts.isSave3D = true;
	nii_loadDir(&opts);
}

void transformDicom2(char *filename, char *output, int processDir)
{
	TDCMopts opts;
	char auxFile[1024], outdir[1024];

	extractFilePath(output, outdir);
	defaultOptions(opts, outdir);

	changeFileExt(output, "", auxFile);
	extractFileName(auxFile, opts.filename);
	if (processDir) 
	{
		extractFilePath(filename, auxFile);
		strcpy(opts.indir, auxFile);
		nii_loadDir(&opts);
	}
	else singleDICOM(&opts, filename);
}