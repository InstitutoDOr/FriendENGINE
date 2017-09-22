#include "dcm2niiInterface.h"
#include "filefuncs.h"

void defaultOptions(TDCMopts &opts, char *outdir)
{
	opts.isCreateBIDS = false;
	opts.isForceStackSameSeries = false;
	opts.isOnlySingleFile = true;
	opts.isCreateText = false;
	opts.isVerbose = 0;
	opts.isCrop = false;
	opts.pigzname[0] = 0;
	opts.filename[0] = 0;
	opts.isGz = false;
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