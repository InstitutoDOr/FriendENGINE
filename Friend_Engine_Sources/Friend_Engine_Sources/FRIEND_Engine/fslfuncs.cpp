#include "newimage/newimageall.h"
#include "libvis/miscplot.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/t2z.h"
#include "libprob.h"
#include "parser.h"
#include "fslio/fslio.h"
#include "intervals.h"
#include "fslfuncs.h"
#include "filefuncs.h"
#include <cstring>
#include <string.h>
#include <vector>

#ifdef WINDOWS
#include <direct.h>
#define strcasecmp stricmp
#define rmdir _rmdir
#endif

using namespace std;
using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

template <class T> int execAxial(char *inNam, char *outNam);
template <class T> int resample(char *inNam, float dx, float dy, float dz, int nn, char *Output);

#ifndef WINDOWS
// functions handling FSL commands
int callCommand(char *command)
{
   char call[10000];
   sprintf(call, "%s/bin/%s", getenv("FSLDIR"), command);
   return system(call);
}

int fslmaths(char *cmdLn)
{
   return callCommand(cmdLn);
}

int fsl_glm(char *cmdLn)
{
   return callCommand(cmdLn);
}

int fslroi(char *cmdLn)
{
   return callCommand(cmdLn);
}

int susan(char *cmdLn)
{
   return callCommand(cmdLn);
}

int bet(char *cmdLn)
{
   return callCommand(cmdLn);
}

int mcflirt(char *cmdLn)
{
   return callCommand(cmdLn);
}

int flirt(char *cmdLn)
{
   return callCommand(cmdLn);
}

int fsl_tsplot(char *cmdLn)
{
   return callCommand(cmdLn);
}

int feat_model(char *cmdLn)
{
   return callCommand(cmdLn);
}

int convert_xfm(char *cmdLn)
{
	return callCommand(cmdLn);
}
#endif

float norm(Matrix &mat)
{
	float normv = 0;
	for (int t = 1; t <= 4; t++)
		normv += mat(t, 1) * mat(t, 1);
	return sqrt(normv);
}

void aff2rigid(char *infile, char *outfile)
{
	Matrix a = read_ascii_matrix(infile);

	// set specific AC and PC coordinates in FLIRT convention(x1 = AC, x2 = PC, x3 = point above x1 in the mid - sag plane)
	Matrix x1(4, 1), x2(4, 1), x3(4, 1);
	x1(1, 1) = 91; x1(2, 1) = 129; x1(3, 1) = 67; x1(4, 1) = 1;
	x2(1, 1) = 91; x2(2, 1) = 100; x2(3, 1) = 70; x2(4, 1) = 1;
	x3(1, 1) = 91; x3(2, 1) = 129; x3(3, 1) = 117; x3(4, 1) = 1;

	Matrix ainv = pinv(a);

	// vectors v are in MNI space, vectors w are in native space
	Matrix v21 = (x2 - x1);
	Matrix v31 = (x3 - x1);

	// normalise and force orthogonality
	v21 = v21 / norm(v21);
	Matrix auxM = v31.t() * v21;
	v31 = v31 - auxM(1, 1) * v21;
	v31 = v31 / norm(v31);

	Matrix v213 = v21.SubMatrix(1, 3, 1, 1);
	Matrix v313 = v31.SubMatrix(1, 3, 1, 1);
	Matrix tmp = cross(v213, v313).t();

	Matrix v41(4, 1);
	v41 = 0;

	for (int t = 1; t <= 3; t++)
		v41(t, 1) = tmp(1, t);

	// Map vectors to native space
	Matrix w21 = ainv*(v21);
	Matrix w31 = ainv*(v31);

	// normalise and force orthogonality
	w21 = w21 / norm(w21);
	auxM = w31.t() * w21;
	w31 = w31 - auxM(1, 1) * w21;
	w31 = w31 / norm(w31);

	v213 = w21.SubMatrix(1, 3, 1, 1);
	v313 = w31.SubMatrix(1, 3, 1, 1);
	tmp = cross(v213, v313).t();

	Matrix w41(4, 1);
	w41 = 0;
	for (int t = 1; t <= 3; t++)
		w41(t, 1) = tmp(1, t);

	// setup matrix : native to MNI space
	Matrix r1(4, 4);
	r1 = 0;
	r1(1, 1) = 1;   	r1(2, 2) = 1;   	r1(3, 3) = 1;   	r1(4, 4) = 1;

	for (int t = 1; t <= 4; t++) r1(t, 1) = w21(t, 1);
	for (int t = 1; t <= 4; t++) r1(t, 2) = w31(t, 1);
	for (int t = 1; t <= 4; t++) r1(t, 3) = w41(t, 1);

	Matrix r2(4, 4);
	r2 = 0;
	r2(1, 1) = 1;   	r2(2, 2) = 1;   	r2(3, 3) = 1;   	r2(4, 4) = 1;

	for (int t = 1; t <= 4; t++) r2(1, t) = v21(t, 1);
	for (int t = 1; t <= 4; t++) r2(2, t) = v31(t, 1);
	for (int t = 1; t <= 4; t++) r2(3, t) = v41(t, 1);

	Matrix r = r2.t() * r1.t();

	// Fix the translation(keep AC = x1 in the same place)
	Matrix ACmni = x1;
	Matrix ACnat = ainv*x1;
	Matrix trans = ACmni - r*ACnat;
	for (int t = 1; t <= 3; t++) r(t, 4) = trans(t, 1);

	// Save out the result
	write_ascii_matrix(r, outfile);
}

void BetFNIRT(char *infile, char *outfile, char *workingDir, char *supportDir)
{
	char CmdLn[2048], fnirtConfig[1024], reference2mm[1024], reference2mmMask[1024], reference[1024], referenceMask[1024], outputBrainMask[1024];
	sprintf(fnirtConfig, "%s/T1_2_MNI152_2mm.cnf", supportDir);

	sprintf(reference2mm, "%s/MNI152_T1_2mm.nii.gz", supportDir);
	sprintf(reference2mmMask, "%s/MNI152_T1_2mm_brain_mask_dil.nii.gz", supportDir);
	sprintf(reference, "%s/MNI152_T1_0.7mm.nii.gz", supportDir);
	sprintf(referenceMask, "%s/MNI152_T1_0.7mm_brain_mask.nii.gz", supportDir);
	sprintf(outputBrainMask, "%s\\bMask.nii", workingDir);

	if (1)
	{
		sprintf(CmdLn, "flirt -interp spline -dof 12 -in %s -ref %s -omat %s/roughlin.mat -out %s/temp_to_MNI_roughlin.nii.gz -nosearch", infile, reference2mm, workingDir, workingDir);
		flirt(CmdLn);
		sprintf(CmdLn, "fnirt --in=%s --ref=%s --aff=%s/roughlin.mat --refmask=%s --fout=%s/str2standard.nii.gz --jout=%s/NonlinearRegJacobians.nii.gz --refout=%s/IntensityModulatedT1.nii.gz --iout=%s/temp_to_MNI_nonlin.nii.gz --logout=%s/NonlinearReg.txt --intout=%s/NonlinearIntensities.nii.gz --cout=%s/NonlinearReg.nii.gz --config=%s", infile, reference2mm, workingDir, reference2mmMask, workingDir, workingDir, workingDir, workingDir, workingDir, workingDir, workingDir, fnirtConfig);
		fnirt(CmdLn);

		// Overwrite the image output from FNIRT with a spline interpolated highres version
		sprintf(CmdLn, "applywarp --rel --interp=spline --in=%s --ref=%s -w %s/str2standard.nii.gz --out=%s/temp_to_MNI_nonlin.nii.gz", infile, reference, workingDir, workingDir);
		applywarp(CmdLn);

		// Invert warp and transform dilated brain mask back into native space, and use it to mask input image
		// Input and reference spaces are the same, using 2mm reference to save time
		sprintf(CmdLn, "invwarp --ref=%s -w %s/str2standard.nii.gz -o %s/standard2str.nii.gz", reference2mm, workingDir, workingDir);
		invwarp(CmdLn);

		sprintf(CmdLn, "applywarp --rel --interp=nn --in=%s --ref=%s -w %s/standard2str.nii.gz -o %s", referenceMask, infile, workingDir, outputBrainMask);
		applywarp(CmdLn);

		sprintf(CmdLn, "fslmaths %s -mas %s %s", infile, outputBrainMask, outfile);
		fslmaths(CmdLn);

		//sprintf(CmdLn, "bet %s %s -f 0.3", outfile, outfile);
		//bet(CmdLn);
	}
	removeDirectory(workingDir);
}

void acpcAlignment(char *inputVolume, char *reference, char *outVolume, char *workingDir, int brainSize, int searchsize)
{
	char CmdLn[2048];
	char infile[1024];
	char oVol[1024];
	char outMatrix[1024];

	sprintf(outMatrix, "%s/saida.mat", workingDir);
	sprintf(oVol, "%s/structural.nii", workingDir);
	strcpy(oVol, inputVolume);
	if (1)
	{
		mkdir(workingDir);
		// Crop the FOV
		sprintf(CmdLn, "robustfov -i %s -m %s/roi2full.mat -r %s/robustroi -b %d", inputVolume, workingDir, workingDir, brainSize);
		robustfov(CmdLn);

		// Invert the matrix(to get full FOV to ROI)
		sprintf(CmdLn, "convert_xfm -omat %s/full2roi.mat -inverse %s/roi2full.mat", workingDir, workingDir);
		convert_xfm(CmdLn);

		// Register cropped image to MNI152(12 DOF)
		// -cost{ mutualinfo, corratio, normcorr, normmi, leastsq, labeldiff, bbr }
		sprintf(CmdLn, "flirt -interp spline -in %s/robustroi.nii.gz -ref %s -omat %s/roi2std.mat -out %s/acpc_final.nii.gz -searchrx -%d %d -searchry -%d %d -searchrz -%d %d", workingDir, reference, workingDir, workingDir, searchsize, searchsize, searchsize, searchsize, searchsize, searchsize);
		flirt(CmdLn);

		// Concatenate matrices to get full FOV to MNI
		sprintf(CmdLn, "convert_xfm -omat %s/full2std.mat -concat %s/roi2std.mat %s/full2roi.mat", workingDir, workingDir, workingDir);
		convert_xfm(CmdLn);

		//Get a 6 DOF approximation which does the ACPC alignment(AC, ACPC line, and hemispheric plane)
		sprintf(infile, "%s/full2std.mat", workingDir);
		aff2rigid(infile, outMatrix);
	}

	if (1)
	{
		// Create a resampled image(ACPC aligned) using spline interpolation
		sprintf(CmdLn, "applywarp --rel --interp=spline -i %s -r %s --premat=%s -o %s", inputVolume, reference, outMatrix, oVol);
		applywarp(CmdLn);
	}

	if (0)
	{
		sprintf(CmdLn, "bet %s %s -f 0.3", oVol, outVolume);
		bet(CmdLn);
		//sprintf(CmdLn, "bet %s %s -f 0.3", outVolume, outVolume);
		//bet(CmdLn);
	}
	removeDirectory(workingDir);
}

void FslFree(FSLIO* OP)
{
	if (OP->niftiptr->fname != NULL) free(OP->niftiptr->fname);
	if (OP->niftiptr->iname != NULL) free(OP->niftiptr->iname);
	if (OP->niftiptr != NULL) free(OP->niftiptr);
	if (OP->fileptr != NULL) free(OP->fileptr);
	free(OP);
}

bool isFSLReadable(char *fileName)
{
	FSLIO *OP = FslOpen(fileName, "r");
	if (OP == NULL) return false;
	else
	{
		FslClose(OP);
		FslFree(OP);
		return true;
	}
}

// maps nifti orientation constants in 1, -1, 2, -2, 3, -3 numbers
int niftiOrientationInt( int ii )
{
   switch( ii ){
	 // X axis
	 case NIFTI_L2R: return 1;
	 case NIFTI_R2L: return -1;
	 // Y axis
	 case NIFTI_P2A: return 2;
	 case NIFTI_A2P: return -2;
	 // Z axis
	 case NIFTI_I2S: return 3;
	 case NIFTI_S2I: return -3;
   }
   return 0;
}

// return ijk orientation from an fslio variable
void returnijk(FSLIO* fslio, int *is, int *js, int *ks)
{
  int i, j, k;
  if (fslio!=NULL)
  {
	 if( fslio->niftiptr->qform_code > 0 ) nifti_mat44_to_orientation( fslio->niftiptr->qto_xyz , &i,&j,&k );
	 else if( fslio->niftiptr->sform_code > 0 ) nifti_mat44_to_orientation( fslio->niftiptr->sto_xyz , &i,&j,&k );
	 else 
	 {
		 i=NIFTI_R2L;
		 j=NIFTI_P2A;
		 k=NIFTI_I2S;

	 }
	 *is = niftiOrientationInt(i);
	*js = niftiOrientationInt(j);
	*ks = niftiOrientationInt(k);
  }
}

// return ijk orientation from a filename
void returnijk(char *inNam, int *is, int *ys, int *zs)
{
  FSLIO* fslio=NULL;

  fslio = FslOpen(FslMakeBaseName(inNam),"rb");
  returnijk(fslio, is, ys, zs);
  if (fslio!=NULL) FslClose(fslio);
}

// get the stats of the clusters
template <class T>
void get_stats(const volume<int>& labelim, const volume<T>& origim,
	vector<int>& size,
	vector<T>& maxvals, vector<float>& meanvals,
	vector<triple<int> >& max, vector<triple<float> >& cog,
	bool minv)
{
	int labelnum = labelim.max();
	size.resize(labelnum + 1, 0);
	maxvals.resize(labelnum + 1, (T)0);
	meanvals.resize(labelnum + 1, 0.0f);
	triple<int> zero;
	zero.x = 0; zero.y = 0; zero.z = 0;
	triple<float> zerof;
	zerof.x = 0; zerof.y = 0; zerof.z = 0;
	max.resize(labelnum + 1, zero);
	cog.resize(labelnum + 1, zerof);
	vector<float> sum(labelnum + 1, 0.0);
	for (int z = labelim.minz(); z <= labelim.maxz(); z++) {
		for (int y = labelim.miny(); y <= labelim.maxy(); y++) {
			for (int x = labelim.minx(); x <= labelim.maxx(); x++) {
				int idx = labelim(x, y, z);
				T oxyz = origim(x, y, z);
				size[idx]++;
				cog[idx].x += ((float)oxyz)*x;
				cog[idx].y += ((float)oxyz)*y;
				cog[idx].z += ((float)oxyz)*z;
				sum[idx] += (float)oxyz;
				if ((size[idx] == 1) ||
					((oxyz>maxvals[idx]) && (!minv)) ||
					((oxyz<maxvals[idx]) && (minv)))
				{
					maxvals[idx] = oxyz;
					max[idx].x = x;
					max[idx].y = y;
					max[idx].z = z;
				}
			}
		}
	}
	for (int n = 0; n <= labelnum; n++) {
		if (size[n]>0.0) {
			meanvals[n] = (sum[n] / ((float)size[n]));
		}
		if (sum[n]>0.0) {
			cog[n].x /= sum[n];
			cog[n].y /= sum[n];
			cog[n].z /= sum[n];
		}
	}
}

void clusterSizeFiltering(char *fileName, char *outputName, int minClusterSize, float minValue, int connectionType)
{
	volume<float> inVol, mask;
	volume<int> componentLabels;
	string inName = string(fileName);
	vector<int> size, idx;
	vector<triple<int> > maxpos;
	vector<triple<float> > cog;
	vector<float> maxvals;
	vector<float> meanvals;

	read_volume(inVol, inName);
	mask = inVol;
	if (minValue != 0) mask.binarise(minValue);
	else mask.binarise(0, mask.max() + 1, exclusive);
	componentLabels = connected_components(mask, connectionType);

	// get the statistics of the connected components
	get_stats<float>(componentLabels, inVol, size, maxvals, meanvals, maxpos, cog, 0);

	// Zeroing the components smaller than minClusterSize
	for (int i = 1; i < size.size(); i++)
	{
		if (size[i] < minClusterSize)
		{
			for (int z = componentLabels.minz(); z <= componentLabels.maxz(); z++){
				for (int y = componentLabels.miny(); y <= componentLabels.maxy(); y++){
					for (int x = componentLabels.minx(); x <= componentLabels.maxx(); x++)
					{
						if (componentLabels(x, y, z) == i) componentLabels(x, y, z) = 0;
					}
				}
			}
		}
	}

	// filtering the volume
	for (int z = componentLabels.minz(); z <= componentLabels.maxz(); z++){
		for (int y = componentLabels.miny(); y <= componentLabels.maxy(); y++){
			for (int x = componentLabels.minx(); x <= componentLabels.maxx(); x++)
			{
				if (componentLabels(x, y, z) == 0) inVol(x, y, z) = 0;
			}
		}
	}
	save_volume(inVol, string(outputName));
}

bool greaterThan(int i, int j) { return (i>j); }
void fillHoles(char *input, char *completeVolume, char *output)
{
	volume<float> inVol, completeVol, mask;
	volume<int> componentLabels;
	vector<int> size, idx, orderedSize;
	vector<triple<int> > maxpos;
	vector<triple<float> > cog;
	vector<float> maxvals;
	vector<float> meanvals;

	read_volume(inVol, string(input));
	read_volume(completeVol, string(completeVolume));
	mask = inVol;
	mask.binarise(0, mask.max() + 1, exclusive);

	for (int z = mask.minz(); z <= mask.maxz(); z++){
		for (int y = mask.miny(); y <= mask.maxy(); y++){
			for (int x = mask.minx(); x <= mask.maxx(); x++)
			{
				if (mask(x, y, z) == 0) mask(x, y, z) = 2;
			}
		}
	}

	componentLabels = connected_components(mask, 26);

	// get the statistics of the connected components
	get_stats<float>(componentLabels, inVol, size, maxvals, meanvals, maxpos, cog, 0);

	if (size.size() > 2)
	{
		for (int i = 0; i < size.size(); i++) orderedSize.push_back(size[i]);
		std::sort(orderedSize.begin(), orderedSize.end(), greaterThan);
		int maxValue = orderedSize[2];
		//fprintf(stderr, "No comps.: %d Corte superior: %d\n", size.size(), maxValue);
		//for (int i = 0; i < size.size(); i++) fprintf(stderr, "%d, ", size[i]);
		//fprintf(stderr, "\n");
		for (int i = 0; i < size.size(); i++)
		{
			if (size[i] <= maxValue)
			{
				for (int z = componentLabels.minz(); z <= componentLabels.maxz(); z++){
					for (int y = componentLabels.miny(); y <= componentLabels.maxy(); y++){
						for (int x = componentLabels.minx(); x <= componentLabels.maxx(); x++)
						{
							if (componentLabels(x, y, z) == i)
								inVol(x, y, z) = completeVol(x, y, z);
						}
					}
				}
			}
		}
	}
	save_volume(inVol, string(output));
}

// template function to reorienting a volume file to axial
template <class T>
int execAxial(char *inNam, char *outNam)
{
  int newx, newy, newz;
  int i, j, k;
  int retVal=0;

  string inName=inNam;
  string outName=outNam;

  volume4D<T> inVol;
  read_orig_volume4D(inVol,inName);
  FSLIO* fslio=NULL;

  fslio = FslOpen(FslMakeBaseName(inNam),"rb");
  if (fslio!=NULL)
  {
	Matrix affMat;
	returnijk(fslio, &i, &j, &k);
	 FslClose(fslio);

	 if (i==1) newx=1;
	 else if (i==-1) newx=-1;
	 else if (j==-1) newx=-2;
	 else if (j==1) newx=2;
	 else if (k==-1) newx=-3;
	 else if (k==1) newx=3;

	 if (i==2) newy=1;
	 else if (i==-2) newy=-1;
	 else if (j==-2) newy=-2;
	 else if (j==2) newy=2;
	 else if (k==-2) newy=-3;
	 else if (k==2) newy=3;

	 if (i==3) newz=1;
	 else if (i==-3) newz=-1;
	 else if (j==-3) newz=-2;
	 else if (j==3) newz=2;
	 else if (k==-3) newz=-3;
	 else if (k==3) newz=3;

	 affMat = inVol.swapmat(newx,newy,newz);
	 float det=affMat.Determinant();
	 for (int t=0; t<inVol.tsize(); t++) {
		inVol[t].swapdimensions(newx,newy,newz);
		if (det < 0) inVol[t].swapdimensions(-1,2,3);
	 }

	 if (outName!="") {
		retVal = save_orig_volume4D(inVol,outName);
	 }
  }
  return retVal;
}

// reorients a volume file to axial
int axial(char *inNam, char *outNam)
{
   // first things first : getting the volume datatype
   string inName=inNam;
   string baseName = fslbasename(inNam);
   
   FSLIO* IP1;
   IP1 = FslOpen(baseName.c_str(),"rb");
   if (IP1==NULL) {
	  cerr << "Cannot open volume " << baseName << " for reading!\n";
   }
   else
   {
	  short type;
	  FslGetDataType(IP1,&type);
	  FslClose(IP1);
	  free(IP1);
	  
	  short dataType=NEWIMAGE::closestTemplatedType(type);
	  
	  // now calling the tempalte function
	  if ( dataType==DT_UNSIGNED_CHAR ) return execAxial<char>(inNam, outNam);
	  else if ( dataType==DT_SIGNED_SHORT ) return execAxial<short>(inNam, outNam);
	  else if ( dataType==DT_SIGNED_INT ) return execAxial<int>(inNam, outNam);
	  else if ( dataType==DT_FLOAT )  return execAxial<float>(inNam, outNam);
	  else if ( dataType==DT_DOUBLE ) return execAxial<double>(inNam, outNam);
   }
   return -1;
}

// applies a threshold in a file based on volume intensities, to create a mask. Not used right now
void brainThreshold(char *functional, char* mask, char *temp, char *output, float brainThres)
{
	float p2, p98;
	float thres;
	char cmd[500];
	char quote[10];
	volume4D<float> Vol;
	
	sprintf(quote, "\"");

   // 2th percentil and 98th percentil calculation
	read_volume4D(Vol, string(functional));
	p2 = Vol.percentile(0.02);
	p98 = Vol.percentile(0.98);
   
   // threshold calculation. Minimum value equals a percentage from the difference between the 98th and 2th percentile. Note the brain threshold parameter is in the range [1..100]
   thres = p2 + (brainThres / 100) * (p98 - p2);
   
   // applying the threshold and get the minimum value per voxel, in case of a more than one volume in file.
   sprintf(cmd, "fslmaths %s%s%s -thr %f -Tmin -bin %s%s%s -odt char", quote, functional, quote, thres, quote, temp, quote);
	fslmaths(cmd);

   // applying a mask in the result
   sprintf(cmd, "fslmaths %s%s%s -mas %s%s%s %s%s%s", quote, mask, quote, quote, temp, quote, quote, output, quote);
	fslmaths(cmd);
}

// returns the intensity value of the perc percentile. Note perc here is in range [0..1]
float percentile(char *volumeFileName, char* maskFileName, float percentage)
{
	volume4D<float> volumeObj;
	volume<float> mask;

	read_volume4D(volumeObj, string(volumeFileName));
	if (maskFileName != NULL)
	{
		read_volume(mask, string(maskFileName));
		return volumeObj.percentile(percentage, mask);
	}
	else return volumeObj.percentile(percentage);
}

// template function to accomodate a volume `mask` to the same  coordinate space of a volume `base`. That is useful in applying a matrix transformation in volumes with different coordinate spaces
template <class T>
int execSameFov(char *base, char *mask, char *output)
{
  short maskOrig[5];
  short baseOrig[5];
  int mx, my, mz;
  Matrix qmat, qbmat;
  volume<T> Base;
  read_volume(Base,string(base));

  volume<T> Mask;
  read_volume(Mask,string(mask));
 

   volume<T>maskSameFov(Base.xsize(),Base.ysize(),Base.zsize());
   maskSameFov.copyproperties(Base);

   // Volume Base center of gravity
   if (Base.sform_code()!=NIFTI_XFORM_UNKNOWN) {
	  qbmat = Base.sform_mat();
	  baseOrig[0]=abs(qbmat(1,4)) + 1;
	  baseOrig[1]=abs(qbmat(2,4)) + 1;
	  baseOrig[2]=abs(qbmat(3,4)) + 1;
   }
   else
   if (Base.qform_code()!=NIFTI_XFORM_UNKNOWN) {
	  qbmat = Base.qform_mat();
	  baseOrig[0]=abs(qbmat(1,4)) + 1;
	  baseOrig[1]=abs(qbmat(2,4)) + 1;
	  baseOrig[2]=abs(qbmat(3,4)) + 1;
   };

   // Volume to be transformed center of gravity
	if (Mask.sform_code()!=NIFTI_XFORM_UNKNOWN) {
	  qmat = Mask.sform_mat();
	  maskOrig[0]=abs(qmat(1,4)) + 1;
	  maskOrig[1]=abs(qmat(2,4)) + 1;
	  maskOrig[2]=abs(qmat(3,4)) + 1;
   }
   else
   if (Mask.qform_code()!=NIFTI_XFORM_UNKNOWN) {
	  qmat = Mask.qform_mat();
	  maskOrig[0]=abs(qmat(1,4)) + 1;
	  maskOrig[1]=abs(qmat(2,4)) + 1;
	  maskOrig[2]=abs(qmat(3,4)) + 1;
   }

   // calculation the transformation mapping between volumes
   Matrix transf = qmat.i() * qbmat;
   float a11=transf(1,1), a12=transf(1,2), a13=transf(1,3), a14=transf(1,4),
		 a21=transf(2,1), a22=transf(2,2), a23=transf(2,3), a24=transf(2,4),
		 a31=transf(3,1), a32=transf(3,2), a33=transf(3,3), a34=transf(3,4);

   // applying the trasformation.
   for(int z=Base.minz(); z<=Base.maxz();z++){
	 for(int y=Base.miny(); y<=Base.maxy();y++){
	   for(int x=Base.minx(); x<=Base.maxx();x++)
	   {
		   mx = x+(maskOrig[0]-baseOrig[0]);
		   my = y+(maskOrig[1]-baseOrig[1]);
		   mz = z+(maskOrig[2]-baseOrig[2]);

		   mx = a11*x + a12*y + a13*z + a14;
		   my = a21*x + a22*y + a23*z + a24;
		   mz = a31*x + a32*y + a33*z + a34;

		   if ((mx >= Mask.minx()) && (mx <= Mask.maxx()) && 
			   (my >= Mask.miny()) && (my <= Mask.maxy()) && 
			   (mz >= Mask.minz()) && (mz < Mask.maxz()))
					maskSameFov(x,y,z) = Mask(mx, my, mz);
		   else maskSameFov(x,y,z) = 0;

	   }
	 }
   }
   // saving the result
   save_volume(maskSameFov, string(output));
   return 1;
}

// accomodate a volume `mask` to the same coordinate space of a volume `base`. Depends on the a template function.
int sameFov(char *base, char *mask, char *output)
{
  string inName=mask;
  string baseName = fslbasename(inName);

  FSLIO* IP1;
  IP1 = FslOpen(baseName.c_str(),"rb");
  if (IP1==NULL) {
	cerr << "Cannot open volume " << baseName << " for reading!\n";
  }
  else
  {
	 // getting the voxel datatype
	 short type;
	 FslGetDataType(IP1,&type);
	 FslClose(IP1);
	 free(IP1);

	 short dataType=NEWIMAGE::closestTemplatedType(type);
	 
	 // calling the template function with the type
	 if ( dataType==DT_UNSIGNED_CHAR ) return execSameFov<char>(base, mask, output);
	 else if ( dataType==DT_SIGNED_SHORT ) return execSameFov<short>(base, mask, output);
	 else if ( dataType==DT_SIGNED_INT ) return execSameFov<int>(base, mask, output);
	 else if ( dataType==DT_FLOAT )  return execSameFov<float>(base, mask, output);
	 else if ( dataType==DT_DOUBLE ) return execSameFov<double>(base, mask, output);
  }
  return -1;
}

// template function to make a `mask` volume in axial and in the same coordinate space and voxel geometry that `base` volume
template <class T>
int execStandardizeVolume(char *base, char *mask, char *output, int NN)
{
	char baseTemp[500];
	volume<T> Base;
	read_volume(Base,string(base));
	sprintf(baseTemp, "%st.nii", FslMakeBaseName(base));
	if (NN == 0)
	{
	  // harmonizes the voxel dimensions
	   resample<float>(mask, Base.xdim(), Base.ydim(), Base.zdim(), NN, output);
	  
	  // axial everybody
	   execAxial<float>(output, output);
	   execAxial<float>(base, baseTemp);
	  
	  // same coordinate space
	   execSameFov<float>(baseTemp, output, output);
	}
	else
	{
	  // harmonizes the voxel dimensions
	   resample<T>(mask, Base.xdim(), Base.ydim(), Base.zdim(), NN, output);
	  
	  // axial everybody
	   execAxial<T>(output, output);
	   execAxial<T>(base, baseTemp);
	  
	  // same coordinate space
	  execSameFov<T>(baseTemp, output, output);
	}

	if (fileExists(baseTemp))
	   remove(baseTemp);
	return 1;
}

// reformats a `mask` volume to axial and harmonizes the voxels geometry with `base` volume. Depends on the above template function.
int standardizeVolume(char *base, char *mask, char *output, int NN)
{
  string inName=mask;
  string baseName = fslbasename(inName);

  FSLIO* IP1;
  IP1 = FslOpen(baseName.c_str(),"rb");
  if (IP1==NULL) {
	cerr << "Cannot open volume " << baseName << " for reading!\n";
  }
  else
  {
	 // getting voxel data type
	 short type;
	 FslGetDataType(IP1,&type);
	 FslClose(IP1);
	 free(IP1);

	 // calling the template function
	 short dataType=NEWIMAGE::closestTemplatedType(type);
	 if ( dataType==DT_UNSIGNED_CHAR ) return execStandardizeVolume<char>(base, mask, output, NN);
	 else if ( dataType==DT_SIGNED_SHORT ) return execStandardizeVolume<short>(base, mask, output, NN);
	 else if ( dataType==DT_SIGNED_INT ) return execStandardizeVolume<int>(base, mask, output, NN);
	 else if ( dataType==DT_FLOAT )  return execStandardizeVolume<float>(base, mask, output, NN);
	 else if ( dataType==DT_DOUBLE ) return execStandardizeVolume<double>(base, mask, output, NN);
  }
  return -1;
}

// calculates the final step of the FRIEND pipeline, the sliding window mean. The `start` variable is just the beginning of the list. For all, always 1. Note this saves a 4D volume
int estimateActivation(int start, int ini, int end, int slidingWindowSize, char *suffix, char *basalFileName, char *maskFileName, char *outputFileName)
{
   char tmpName[5000];

   // verifying if all files exists
   for (int i = ini; i <= end; i++)
   {
	   sprintf(tmpName, suffix, ini);
	   if (!fileExists(tmpName)) // file Exists corrects the name, like adding .gz, etc..
	   {
		   fprintf(stderr, "File %s does not exists. Leaving the estimateActivation function.\n", tmpName);
		   return 0;
	   }
   }

   // reading basal volume, if exists
   volume<float> basalVolume;
   if (basalFileName != NULL)
   {
	  string bName = basalFileName;
	  read_volume(basalVolume, bName);
   }

   // suffix here handles all the the width of the number, e.g 00001, not just 1. It cames like %.5d
   volume<float> auxVol;
   sprintf(tmpName, suffix, ini);

   // until the file is ready, stall. That's the best option?
   if (fileExists(tmpName)) // file Exists corrects the name, like adding .gz, etc..
	  while (isFSLReadable(tmpName)==false) {};
   string auxn = tmpName;
   
   
   // read first volume just to get the volume properties
   read_volume(auxVol, auxn);
   volume4D<float>activation(auxVol.xsize(),auxVol.ysize(),auxVol.zsize(),end-ini+1);
   
   // and copy to the output
   activation.copyproperties(auxVol);
   
   // iterates for each volume
   for(int i=ini; i <= end; i++)
   {
	   int beginning=Max(start, i-slidingWindowSize+1);
	   activation[i-ini] = 0;
	  // sliding window mean
	   for (int j=beginning; j<=i; j++)
	   {
		  volume<float> tmp;
		  sprintf(tmpName, suffix, j);
		  if (fileExists(tmpName)) // file Exists corrects the name, like adding .gz, etc.. 
			 while (isFSLReadable(tmpName)==false) {};
		  string Tmpname = tmpName;
		  read_volume(tmp, Tmpname);
		  activation[i-ini] += tmp;
	   }
	   activation[i-ini] /= (float) (i-beginning+1);
	  // if `basalFileName` is defined, then subtract the value
	   if (basalFileName != NULL) activation[i-ini]-=basalVolume;
   }
   
   // if maskFileName is defined, then apply then mask in the output
   if (maskFileName != NULL)
   {
	   volume<float> mask;
	   string maskFile = maskFileName;
	   read_volume(mask, maskFile);
	   mask.binarise(0,mask.max()+1,exclusive);
	   for (int i=ini; i<=end; i++) activation[i-ini]*=mask;
   }
   
   // saving the result
   string outName = outputFileName;
   save_volume4D(activation, outName); 
   return 0;
}

// same function above, but with default values for start and basalFilename
int estimateActivation(int ini, int end, int slidingWindowSize, char *suffix, char *maskFileName, char *output)
{
	return estimateActivation(1, ini, end, slidingWindowSize, suffix, NULL, maskFileName, output);
}

// same function above, but with default values for start, basalFilename and maskFilename
int estimateActivation(int ini, int end, int slidingWindowSize, char *suffix, char *output)
{
   return estimateActivation(1, ini, end, slidingWindowSize, suffix, NULL, NULL, output);
}

// the following in memory functions are not used right now due to no existing implementations of mcflirt and flirt like in Windows

// loads a mask to memory
int memOpenBinarizedMask(char *Mask, volume<float> *&mask)
{
	mask = new volume<float>;
	string maskFile = Mask;
	read_volume(*mask, maskFile);
	mask->binarise(0,mask->max()+1,exclusive);
	return 0;
}

// calculates the final step of the FRIEND pipeline, the sliding window mean with all data in memory for just one volume.
int memEstimate4DVolumeActivation(int start, int ini, int slidingWindowSize, volume4D<float> *vol4D, volume<float> *basal, volume<float> *mask, volume<float> *&activation)
{
	activation = new volume<float>;
	activation->reinitialize(vol4D->xsize(),vol4D->ysize(),vol4D->zsize());
	activation->copyproperties((*vol4D)[0]);
	int beginning=Max(start, ini-slidingWindowSize+1);
	*activation = 0;
	for (int j=beginning; j<=ini; j++) *activation += (*vol4D)[j-start];
	*activation /= (float) (ini-beginning+1);
	if (basal != NULL) *activation-=*basal;

	if (mask != NULL) (*activation)*=(*mask);
	return 0;
}

// calculates the final step of the FRIEND pipeline, the sliding window mean with all data in memory, except for basal and mask volumes
int memEstimate4DActivation(int start, int ini, int end, int slidingWindowSize, volume4D<float> *vol4D, char *basalFileName, char *maskFileName, char *output)
{
   int inc=start;

   // loads the basal file, if defined
   volume<float> basalVolume;
   if (basalFileName != NULL)
   {
	  string bName = basalFileName;
	  read_volume(basalVolume, bName);
   }
   
   volume4D<float>activation(vol4D->xsize(),vol4D->ysize(),vol4D->zsize(),end-ini+1);
   activation.copyproperties(*vol4D);
   
   // calculates the sliding window mean
   for(int i=ini; i <= end; i++)
   {
	   int beginning=Max(start, i-slidingWindowSize+1);
	   activation[i-ini] = 0;
	   for (int j=beginning; j<=i; j++) activation[i-ini] += (*vol4D)[j-inc];
	   activation[i-ini] /= (float) (i-beginning+1);
	   if (basalFileName != NULL) activation[i-ini]-=basalVolume;
   }

   // apply a mask volume in the result, if defined
   if (maskFileName != NULL)
   {
	   volume<float> maskVolume;
	   string maskFile = maskFileName;
	   read_volume(maskVolume, maskFile);
	   maskVolume.binarise(0,maskVolume.max()+1,exclusive);
	  activation *=maskVolume;
   }

   // save the results
   string outName = output;
   save_volume4D(activation, outName); 
   return 0;
}

// like above function, except that the we have only the 4D volume filename. Here we just open it.
int estimateActivation4D(int start, int ini, int end, int slidingWindowSize, char *vol4DFileName, char *basalFileName, char *maskFileName, char *output)
{
   volume4D<float> vol4D;
   if (vol4DFileName != NULL)
   {
	  string file4DVol = vol4DFileName;
	  read_volume4D(vol4D, file4DVol);
   }
   return memEstimate4DActivation(start, ini, end, slidingWindowSize, &vol4D, basalFileName, maskFileName, output);
}

// the main difference from this and fslroi command is that here we have a list, like, (1, 7, 7, 13, 9) not just a range of indexes.
int filter4D(char *vol4DFileName, vector <int> &indexes, char *output)
{
   // opening the 4d volume
   volume4D<float> vol4D;
   if (vol4DFileName != NULL)
   {
	  string file4DVol = vol4DFileName;
	  read_volume4D(vol4D, file4DVol);
   }

   volume4D<float>activation(vol4D.xsize(),vol4D.ysize(),vol4D.zsize(), indexes.size());
   activation.copyproperties(vol4D);
   
   // picking the indexes we want
   for(int i=0; i < indexes.size(); i++)  activation[i] = vol4D[indexes[i]-1];
   
   // saving the 4d volume
   string outName = output;
   save_volume4D(activation, outName); 
   return 0;
}

// extract the rotation and translations params from a matrix
int returnParams(char *matrixFile, char *referenceFileName, float *rx, float *ry, float *rz, float *tx, float *ty, float *tz)
{
   Matrix affMat(4,4);
   ColumnVector params(12), cor(3);
   cor=0;  // centre of rotations
   
   // reads the matrix
   affMat = read_ascii_matrix(matrixFile);
   if (affMat.Nrows()<4) return -2;
   
   // reads the centre of rotations from the reference volume, if defined
   if (referenceFileName)
   {
	  string volName=referenceFileName;
	  volume<float> reference;
	  if (read_volume(reference,volName)<0)  return -1;
	  cor = reference.cog("scaled_mm");
   };
   // get the parameters
   decompose_aff(params,affMat,cor,rotmat2euler);

   // returning them
   *rx=params(1);
   *ry=params(2);
   *rz=params(3);

   *tx=params(4);
   *ty=params(5);
   *tz=params(6);
   return 0;
}

// extract the rotation and translations params and rms values from a matrix and writes the result in two separeted files.
void matrixParameters(char *matrixFileName, char *referenceFileName, char *output)
{
   float rmax = 80.0;
   ColumnVector vecParam(12);
   ColumnVector center(3);
   Matrix matrix = read_ascii_matrix(string(matrixFileName));
   
   volume<float>reference;
   read_volume(reference, string(referenceFileName));
   
   ofstream outFile, rmsAbsFile;
   string fileName = string(output) + ".par";
   string rmsAbsFileName = string(output) + "_abs.rms";
   
   vecParam = 0;
   center(1) = 0.5*(reference.xsize() - 1.0)*reference.xdim();
   center(2) = 0.5*(reference.ysize() - 1.0)*reference.ydim();
   center(3) = 0.5*(reference.zsize() - 1.0)*reference.zdim();
   
   rmsAbsFile. open(rmsAbsFileName.c_str());
   rmsAbsFile << rms_deviation(IdentityMatrix(4), matrix, center, rmax) << endl;
   
   outFile. open(fileName.c_str());
   decompose_aff(vecParam, matrix, reference.cog("scaled_mm"), rotmat2euler);
   outFile << vecParam(1) << "  " << vecParam(2) << "  "
   << vecParam(3) << "  " << vecParam(4) << "  "
   << vecParam(5) << "  " << vecParam(6) << "  " << endl;
}

// loads the data struture of a volume in memory
FSLIO * fslioopen(char *file)
{
  FSLIO *src = FslOpen(file, "r");
  return src;  
}

// frees the allocated memory
void fslioclose(FSLIO *src)
{
  FslClose(src);
}

// handles volume file axis inversions, if needed, like analyze files from Philips DRINDumper
template <class T>
int execFslSwapDimRT(int argc,char *argv[], FSLIO *src)
{
  string newx=argv[2], newy=argv[3], newz=argv[4];
  string inName=argv[1];
  string outName="";

  if (argc>5) {
	outName=argv[5];
  };

  volume4D<T> inVol;
  read_orig_volume4D(inVol,inName);

  Matrix affMat;
  affMat = inVol.swapmat(newx,newy,newz);
  
  inVol.swapdimensions(newx,newy,newz);

  int retVal=0;

  if (src!=NULL) {
	  float vx, vy, vz, tr;
	  short scode, qcode;
	  mat44 smat, qmat;
	  Matrix sMat(4, 4), qMat(4,4);

	  {
		 scode = FslGetStdXform(src,&smat);
		 qcode = FslGetRigidXform(src,&qmat);
		 FslGetVoxDim(src, &vx, &vy, &vz, &tr);
		 for (int i=0; i<4; i++)
		 {
			 for (int j=0; j<4;j++)
			 {
				sMat.element(i,j) = smat.m[i][j];
				qMat.element(i,j) = qmat.m[i][j];
			 }
		 }

		 inVol.set_qform(qcode, qMat);
		 inVol.set_sform(scode, sMat);
		 inVol.setxdim(vx);
		 inVol.setydim(vy);
		 inVol.setzdim(vz);
		 inVol.setTR(tr);
	  }
  }

  if (outName!="") {
	retVal = save_orig_volume4D(inVol,outName);
  }
  return retVal;
}

// same functionality of the fslswapdim command
int fslSwapDimRT(const char *cmdLn, FSLIO *src)
{
   int r;
   int argc;
   char **argv;
   
   parser(cmdLn, argc, argv);
   
   string progName=argv[0];
   if (argc<5) {
	  //print_usage(progName);
	  r= -1;
   }
   else
   {
	  string inName = argv[1];
	  short dataType1=NEWIMAGE::closestTemplatedType(dtype(inName));
	  if ( dataType1==DT_UNSIGNED_CHAR ) r=execFslSwapDimRT<char>(argc, argv, src);
	  else if ( dataType1==DT_SIGNED_SHORT ) r=execFslSwapDimRT<short>(argc, argv, src);
	  else if ( dataType1==DT_SIGNED_INT ) r=execFslSwapDimRT<int>(argc, argv, src);
	  else if ( dataType1==DT_FLOAT )  r=execFslSwapDimRT<float>(argc, argv, src);
	  else if ( dataType1==DT_DOUBLE ) r=execFslSwapDimRT<double>(argc, argv, src);
	  else r=-1;
   }
   
   freeparser(argc, argv);
   return r;
}

// same as above, but with as variables in memory, except the input volume
int memFslSwapDim_rt(char *cmdLn, FSLIO *src, volume4D<float>* &output)
{
   int argc;
   char **argv;
  
   parser(cmdLn, argc, argv);

   string newx=argv[2], newy=argv[3], newz=argv[4];
   string inName=argv[1];

   output = new volume4D<float>;
   read_orig_volume4D(*output,inName);

   Matrix affMat;
   affMat = output->swapmat(newx,newy,newz);

   output->swapdimensions(newx,newy,newz);
   if (affMat.Determinant() < 0.0)
	  output->swapdimensions(-1,2,3);

   int retVal=0;

   if (src!=NULL) 
   {
	  float vx, vy, vz, tr;
	  short scode, qcode;
	  mat44 smat, qmat;
	  Matrix sMat(4, 4), qMat(4,4);
	  scode = FslGetStdXform(src,&smat);
	  qcode = FslGetRigidXform(src,&qmat);
	  FslGetVoxDim(src, &vx, &vy, &vz, &tr);
	  for (int i=0; i<4; i++)
	  {
		 for (int j=0; j<4;j++)
		 {
			sMat.element(i,j) = smat.m[i][j];
			qMat.element(i,j) = qmat.m[i][j];
		 }
	  }
	  output->set_qform(qcode, qMat);
	  output->set_sform(scode, sMat);
	  output->setxdim(vx);
	  output->setydim(vy);
	  output->setzdim(vz);
	  output->setTR(tr);
   }
   freeparser(argc, argv);
   return retVal;
}


// template function change voxel dimensions of a volume
template <class T>
  int resample(char *inNam, float dx, float dy, float dz, int nn, char *output)
  {
	string inName=inNam;

	volume<T> aniso;
	read_orig_volume(aniso,inName);

	extrapolation oldex = aniso.getextrapolationmethod();
	if ((oldex==boundsassert) || (oldex==boundsexception)) 
	{ 
		aniso.setextrapolationmethod(constpad);
	}

	float stepx, stepy, stepz;
	if (nn==1) aniso.setinterpolationmethod(nearestneighbour);
	else aniso.setinterpolationmethod(trilinear);

	stepx = dx / aniso.xdim();
	stepy = dy / aniso.ydim();
	stepz = dz / aniso.zdim();
	int sx, sy, sz;
	sz = (int) Max(1.0, ( ((float) (aniso.maxz() - aniso.minz() + 1.0)) / stepz));
	sy = (int) Max(1.0, ( ((float) (aniso.maxy() - aniso.miny() + 1.0)) / stepy));
	sx = (int) Max(1.0, ( ((float) (aniso.maxx() - aniso.minx() + 1.0)) / stepx));
	volume<T> iso(sx,sy,sz);
	float fx, fy, fz;
	int x, y, z;
	for (fz=0.0, z=0; z<sz; z++, fz+=stepz) 
	{
		for (fy=0.0, y=0; y<sy; y++, fy+=stepy) 
		{
			for (fx=0.0, x=0; x<sx; x++, fx+=stepx) 
			{
			   iso(x,y,z) = (T)aniso.interpolate(fx,fy,fz);
			}
		}
	}
	iso.copyproperties(aniso);
	iso.setdims(dx,dy,dz);
	 
	// transform the sform and qform matrix appropriately (if set)
	Matrix iso2Aniso(4,4);
	iso2Aniso = 0.0;
	iso2Aniso(1,1)=stepx;
	iso2Aniso(2,2)=stepy;
	iso2Aniso(3,3)=stepz;
	iso2Aniso(4,4)=1.0;

	if (aniso.sform_code()!=NIFTI_XFORM_UNKNOWN) {
		Matrix result = aniso.sform_mat() * iso2Aniso;
	  iso.set_sform(aniso.sform_code(), result);
	}
	if (aniso.qform_code()!=NIFTI_XFORM_UNKNOWN) {
		Matrix result = aniso.qform_mat() * iso2Aniso;
		iso.set_qform(aniso.qform_code(), result);
	}
	string fileName=output;
	save_orig_volume(iso, fileName);
	aniso.setextrapolationmethod(oldex);
	return 0;
}

// changes the voxel dimension of a volume. Depends of resample function
int resampleVolume(char *anatomic, char *sAnatomic, float dx, float dy, float dz, float TR, int nn)
{
   FSLIO* IP1;
   IP1 = FslOpen(anatomic,"rb");
   if (IP1==NULL) {
	  cerr << "Cannot open volume " << anatomic << " for reading!\n";
   }
   else
   {
	  short type;
	  FslGetDataType(IP1,&type);
	  FslClose(IP1);
	  free(IP1);
	  
	  short dataType=NEWIMAGE::closestTemplatedType(type);
	  if (nn ==0) return resample<float>(anatomic, dx, dy, dz, nn, sAnatomic);
	  else
	  {
		 if ( dataType==DT_UNSIGNED_CHAR ) return resample<char>(anatomic, dx, dy, dz, nn, sAnatomic);
		 else if ( dataType==DT_SIGNED_SHORT ) return resample<short>(anatomic, dx, dy, dz, nn, sAnatomic);
		 else if ( dataType==DT_SIGNED_INT ) return resample<int>(anatomic, dx, dy, dz, nn, sAnatomic);
		 else if ( dataType==DT_FLOAT )  return resample<float>(anatomic, dx, dy, dz, nn, sAnatomic);
		 else if ( dataType==DT_DOUBLE ) return resample<double>(anatomic, dx, dy, dz, nn, sAnatomic);
	  }
   }
   return -1;
   
}

// resample `volume` file to the same voxel dimensions that `reference`
void equalVoxelDim(char *volumeFileName, char *reference, char *output, float TR, int nn)
{
   volume<float> ref;
   string refFile=reference;
   read_orig_volume(ref, refFile);
   resampleVolume(volumeFileName, output, ref.xdim(), ref.ydim(), ref.zdim(), TR, nn);
}

// template function to centralize volume in FOV
template <class T>
int execCentralizeVolume(char *volumeFileName, char *output)
{
	string inName=volumeFileName;
	volume<T> aniso;
	read_orig_volume(aniso,inName);
	int dimMax = Max(aniso.zsize(), Max(aniso.ysize(), aniso.xsize()));
	volume<T> newVolume(dimMax, dimMax, dimMax);
	pad(aniso, newVolume);

	string fileName=output;
	save_orig_volume(newVolume, fileName);
	return 1;
}

// centralizes a volume in FOV
int centralizeVolume(char *volume, char *sAnatomic)
{
  FSLIO* IP1;
  IP1 = FslOpen(volume,"rb");
  if (IP1==NULL) {
	cerr << "Cannot open volume " << volume << " for reading!\n";
  }
  else
  {
	 short type;
	 FslGetDataType(IP1,&type);
	 FslClose(IP1);
	 free(IP1);

	 short dataType=NEWIMAGE::closestTemplatedType(type);
	 if ( dataType==DT_UNSIGNED_CHAR ) return execCentralizeVolume<char>(volume, sAnatomic);
	 else if ( dataType==DT_SIGNED_SHORT ) return execCentralizeVolume<short>(volume, sAnatomic);
	 else if ( dataType==DT_SIGNED_INT ) return execCentralizeVolume<int>(volume, sAnatomic);
	 else if ( dataType==DT_FLOAT )  return execCentralizeVolume<float>(volume, sAnatomic);
	 else if ( dataType==DT_DOUBLE ) return execCentralizeVolume<double>(volume, sAnatomic);
  }
  return -1;
}

// transforms `mniTemplate` volume in MNI space to native space of `betRFI`
void MniToSubject(char *betRFI, char * mniTemplate, char * mniStandard, char* RFI2MNI, char* RFI2MNITransf, char* MNI2RFITransf, char* output, char *prefix)
{
	char interpolation[100];
	char standardTemplate[500];
	char prefixcalc[500];
	char outDir[500];
	char ext[]="_std.nii";
	stringstream cmdLn;
	
	strcpy(interpolation, " -interp nearestneighbour");
	if (prefix==NULL) strcpy(prefixcalc, "_RFI2");
	else strcpy(prefixcalc, prefix);
   
	// make standard and template in the same bounding box
	extractFilePath(betRFI, outDir);
	includeTrailingPathDelimiter(outDir);
	changeFileExt(mniTemplate, ext, standardTemplate);
	if (!fileExists(standardTemplate))
		standardizeVolume(mniStandard, mniTemplate, standardTemplate, 1);
	
	// co-register RFI with standard
	sprintf(RFI2MNITransf, "%s%s%s%s", outDir, "RFI2MNI", prefixcalc, ".mat");
	sprintf(RFI2MNI, "%s%s%s%s", outDir, "RFI2MNI", prefixcalc, ".nii");
	if (!fileExists(RFI2MNITransf))
	{
		cmdLn.str("");
		cmdLn << "flirt -in \"" << betRFI << "\" -ref \"" << mniStandard << "\" -o \"" << RFI2MNI << "\" -omat \"" << RFI2MNITransf << "\"";
		flirt((char *)cmdLn.str().c_str());
	};
	
	// calculate the inverse transformation
	sprintf(MNI2RFITransf, "%s%s%s%s", outDir, "MNI2RFI", prefixcalc, ".mat");
	if (!fileExists(MNI2RFITransf))
	{
		cmdLn.str("");
		cmdLn << "convert_xfm -inverse -omat \"" << MNI2RFITransf << "\" \"" << RFI2MNITransf << "\"";
		convert_xfm((char *)cmdLn.str().c_str());
	};
	
	// apply inverse into Template to bring to RFI space
	if (!fileExists(output))
	{
		cmdLn.str("");
		cmdLn << "flirt -ref \"" <<  betRFI << "\" -in \"" << standardTemplate << "\" -o \"" << output << "\" -applyxfm -init \"" << MNI2RFITransf << "\" " << interpolation;
		flirt((char *)cmdLn.str().c_str());

		// make sure of the right side. Further calculations will use the volume without using the transformation matrix
		sameFov(betRFI, output, output);
	};
};

// transforms `mniTemplate` volume in MNI space to native space of `betRFI`
void MniToSubject(char *betRFI, char *mniTemplate, char *mniStandard, char* output, char *prefix)
{
	char RFI2MNI[500], RFI2MNITransf[500], MNI2RFITransf[500];
   
	MniToSubject(betRFI, mniTemplate, mniStandard, RFI2MNI, RFI2MNITransf, MNI2RFITransf, output, prefix);
}

// this function adjusts a `mask` volume corregistered with a `functional` with another `reference` volume. This function is used to adjust a mask coregistered with the RFI volume with the first volume of a run, to account for movements.
void functionalNormalization(char *mask, char *functional, char *reference, char *output, bool nearestNeighbour)
{
   char matOldName[BUFF_SIZE], matNewName[BUFF_SIZE];
   if (strcasecmp(functional, reference) == 0) copyfile(mask, output);
   else
   {
	  stringstream cmdLn;
	  cmdLn << "mcflirt -in " << functional << " -reffile " << reference << " -out " << output << " -mats";
	  mcflirt((char *)cmdLn.str().c_str());
	  
	  // copying the .mat file
	  // from file
	  sprintf(matOldName, "%s%s%c%s", output, ".mat", PATHSEPCHAR, "MAT_0000");
	  // to file
	  changeFileExt(output, ".txt", matNewName);
	  rename(matOldName, matNewName);

	  // deleting directory
	  sprintf(matOldName, "%s%s", output, ".mat");
	  rmdir(matOldName);

	  // renaming from .txt to .mat. WE only can do that after deleting the directory
	  changeFileExt(output, ".txt", matOldName);
	  changeFileExt(output, ".mat", matNewName);
	  rename(matOldName, matNewName);

	  
	  // for now the nearest neighbour is set true no matter what
	  cmdLn.str("");
	  cmdLn << "flirt -in " << mask << " -ref " << reference << " -out " << output << " -applyxfm -paddingsize 10 -init " << matNewName <<
	  " -interp nearestneighbour";
	  flirt((char *)cmdLn.str().c_str());
   }
}

// this function engraves a roi volume in a RFI volume, to make sure of the side
void uniteVolumes(char *referenceVolume, char *roiVolume, char *outputFile)
{
	volume<float> reference, roi;

	read_volume(reference, string(referenceVolume));
	read_volume(roi, string(roiVolume));

	for (int z = roi.minz(); z <= roi.maxz(); z++){
		for (int y = roi.miny(); y <= roi.maxy(); y++){
			for (int x = roi.minx(); x <= roi.maxx(); x++)
			{
				if (roi(x, y, z)) reference(x, y, z) = roi(x, y, z);
			}
		}
	}
	save_volume(reference, string(outputFile));
}

void invertIndexes(vector<int> &idxs, int size)
{
	vector<int>temp;
	temp = idxs;
	idxs.clear();
	std::sort(temp.begin(), temp.end());

	int idxVector = 0;
	for (int t = 1; t <= size; t++)
	{
		if (idxVector < temp.size())
		{
			if (t < temp[idxVector]) idxs.push_back(t);
			else idxVector++;
		}
		else idxs.push_back(t);
	}
}

void generateTMaxVoxels(char *input, char *temp, char *output, vector<int> &idxs, int invert)
{
	volume4D<float>inputVolume;
	volume4D<float>outputVolume;

	read_volume4D(inputVolume, string(input));
	if (invert) invertIndexes(idxs, inputVolume.tsize());

	outputVolume.reinitialize(inputVolume.xsize(), inputVolume.ysize(), inputVolume.zsize(), idxs.size());
	outputVolume.copyproperties(inputVolume);

	for (int t = 0; t < idxs.size(); t++)
	{
		outputVolume[t] = inputVolume[idxs[t] - 1];
	}

	save_volume4D(outputVolume, string(temp));
	stringstream CmdLn;
	CmdLn << "fslmaths " << temp << " -Tmax " << output;
	fslmaths((char *)CmdLn.str().c_str());
}

void generateTMaxVoxels(char *input, char *output, vector<int> &idxs, int invert)
{
	volume4D<float>inputVolume;
	volume<float>tmaxVolume;

	read_volume4D(inputVolume, string(input));
	if (invert) invertIndexes(idxs, inputVolume.tsize());

	tmaxVolume.reinitialize(inputVolume[0]);
	tmaxVolume.copyproperties(inputVolume[0]);
	int numVoxels = tmaxVolume.xsize() * tmaxVolume.ysize() * tmaxVolume.zsize();

	float *tmaxDataPtr = (float *)tmaxVolume.fbegin();

	for (int t = 0; t < numVoxels; t++)
	{
		float maxValue = -10000;
		for (int j = 0; j < idxs.size(); j++)
		{
			maxValue = MAX(maxValue, inputVolume[idxs[j] - 1].fbegin()[t]);
		}
		tmaxDataPtr[t] = maxValue;
	}
	save_volume(tmaxVolume, string(output));
}
