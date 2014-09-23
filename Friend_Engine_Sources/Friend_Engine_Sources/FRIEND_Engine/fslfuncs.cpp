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

using namespace std;
using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

#ifdef WINDOWS
#include <direct.h>
#define strcasecmp stricmp
#define rmdir _rmdir
#endif

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
        //if (det < 0) inVol[t].swapdimensions(-1,2,3);
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
	sprintf(baseTemp, "%st.nii.gz", FslMakeBaseName(base));
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
   string auxn = tmpName;

   // until the file is ready, stall. That's the best option?
   while (isReadable(tmpName)==false) {};
   
   
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
         while (isReadable(tmpName)==false) {};
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
	write_ascii_matrix(iso2Aniso, "e:\\matrixiso.txt");

    if (aniso.sform_code()!=NIFTI_XFORM_UNKNOWN) {
		write_ascii_matrix(aniso.sform_mat(), "e:\\matrixaniso.txt");
		Matrix result = aniso.sform_mat() * iso2Aniso;
		write_ascii_matrix(result, "e:\\matrixr1.txt");
      iso.set_sform(aniso.sform_code(), result);
    }
    if (aniso.qform_code()!=NIFTI_XFORM_UNKNOWN) {
		write_ascii_matrix(aniso.qform_mat(), "e:\\matrixaniso.txt");
		Matrix result = aniso.qform_mat() * iso2Aniso;
		write_ascii_matrix(result, "e:\\matrixr2.txt");
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
    char ext[]="_std.nii.gz";
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
      changeFileExt(output, ".mat", matNewName);
      sprintf(matOldName, "%s%s%c%s", output, ".mat", PATHSEPCHAR, "MAT_0000");
      rename(matOldName, matNewName);
      sprintf(matOldName, "%s%s", output, ".mat");
      rmdir(matOldName);

      
      cmdLn.str("");
      cmdLn << "flirt -in " << mask << " -ref " << reference << " -out " << output << " -applyxfm -paddingsize 10 -init " << matNewName <<
      " -interp nearestneighbour";
      flirt((char *)cmdLn.str().c_str());
   }
}