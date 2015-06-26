#include "svm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "svmfuncs.h"
#include "fslfuncs.h"
#include "filefuncs.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/fsl_isfinite.h"
#include "libprob.h"
#include "parser.h"
#include <string>

using namespace MISCMATHS;
using namespace NEWIMAGE;


#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void generateProjetionsGraph(char *projectionsFile)
{
	char pngFile[1000];
	changeFileExt(projectionsFile, ".png", pngFile);
	generateProjetionsGraph(projectionsFile, pngFile);
}

void generateProjetionsGraph(char *projectionsFile, char *pngFile)
{
	char tempFile[1000];
	stringstream CmdLn;
	string row;

	fstream projections(projectionsFile, fstream::in);
	changeFileExt(projectionsFile, ".proj", tempFile);
	fstream Output(tempFile, fstream::in | fstream::out | fstream::trunc);
	while (getline(projections, row))
	{
		std::size_t found = row.find("Accuracy");
		if (found == std::string::npos)
			Output << row << '\n';
	}
	Output.close();
	projections.close();

	// generating the projections graph png
	fprintf(stderr, "Generating svm projection graphics\n");
	CmdLn << "fsl_tsplot -i " << tempFile << " -t \"SVM projections\" -u 1 --start=2 --finish=2 -a Projections -w 640 -h 144 -o " << pngFile;

	fsl_tsplot((char *)CmdLn.str().c_str());
}

// calculates the svm weights of a model, storing the values in a vector object
void calculateWeightVector(const svm_model *model, vector <double> &weightVector)
{
	int l = model->l;
	int nr_class = model->nr_class;
	double coef;
	const double * const *sv_coef = model->sv_coef;
	const svm_node * const *SV = model->SV;

	for (int i=0; i < weightVector.size(); i++) weightVector[i] = 0;

	for(int i=0; i<l; i++)
	{
		coef = 0;
		for(int j=0;j < nr_class-1;j++)
			coef += sv_coef[j][i];

		const svm_node *p = SV[i];
		while ((p->index != -1) && (p->index <= weightVector.size()))
		{
			weightVector[p->index-1] += coef * p->value;
			p++;
		}
	}
}

// get a svm score based on a vector sample. This functions returns the prediction class as a return value
double svmSamplePredict(const svm_model *model, const svm_node *sample, double &score)
{
   int nr_class = model->nr_class;
   double *dec_values;
   if(model->param.svm_type == ONE_CLASS ||
	   model->param.svm_type == EPSILON_SVR ||
	   model->param.svm_type == NU_SVR)
   dec_values = Malloc(double, 1);
   else dec_values = Malloc(double, nr_class*(nr_class-1)/2);
   
   double classPrediction = svm_predict_values(model, sample, dec_values);
   score = dec_values[0];
   free(dec_values);
   return classPrediction;
}

// loads a model file in memory
svm_model* loadModel(char *modelFileName)
{
   svm_model *model;
   if((model=svm_load_model(modelFileName))==0) return NULL;
   else return model;
}

// unloads a model in memory
void unloadModel(svm_model *&model)
{
   svm_free_and_destroy_model(&model);
   model = NULL;
}

// returns the number of voxels whose intensities are above `minValue`
int getVolumeSize(const char *volumeFileName, float minValue)
{
	volume<float> volumeObject;
	int volumeSize=0;
    if (volumeFileName != NULL)
    {
      string volumeFile = volumeFileName;
      read_volume<float>(volumeObject, volumeFile);
      for(int z=0; z < volumeObject.zsize(); z++)
         for(int y=0; y < volumeObject.ysize(); y++)
	        for(int x=0; x < volumeObject.xsize(); x++)
                if (volumeObject.value(x,y,z) > minValue) volumeSize++;
    }
	return volumeSize;
}

// transforms a 4D Volume in a SVM samples file, based on a mask
void saveSVMFile(const char *volume4DFileName, const char *maskFileName, const char *outputFileName, float minValue, vector <int > &indexes, vector <int> &classes)
{
   // reading the mask file
	volume<float> mask;
   if (maskFileName != NULL)
   {
      string Maskfile = maskFileName;
      read_volume<float>(mask, Maskfile);
   }

   // reading the 4D volume
	volume4D<float> volSamples;
   if (volume4DFileName != NULL)
   {
      string Vol4DFile = volume4DFileName;
      read_volume4D<float>(volSamples, Vol4DFile);
   }

	FILE *f;
	f=fopen(outputFileName, "wt+");
	if (f !=NULL)
	{
	   int i, t;
      for(int h=0;h < indexes.size();h++)
	   {
        // picking the right indexes
		  t = indexes[h] - 1;
		  i = 0;
        // saving the class value. Remember indexes size is different from classes. classes has the same sime of the number of volumes
		  fprintf(f, "%d ", classes[t]);
		  for(int z=0; z < mask.zsize(); z++)
			 for(int y=0; y < mask.ysize(); y++)
				for(int x=0; x < mask.xsize(); x++)
                if (mask.value(x,y,z) > minValue)
			       {
					   i++;
                  // writing each voxel value in svm format
					   fprintf(f, "%d:%f ", i, volSamples.value(x,y,z,t));
			       }
		  fprintf(f, "\n");
	   }
	   fclose(f);
	}
}

// functions that transforms a volume into a svm sample vector
void _volume2Sample(svm_model *model, volume<float> &volSample, volume<float> &mask, int sampleSize, float minValue, svm_node * &sample)
{	
	sample=(struct svm_node *) malloc((sampleSize+1)*sizeof(struct svm_node));
    
	int i=0;
   for(int z=0; z < volSample.zsize();z++)
      for(int y=0; y < volSample.ysize();y++)
	     for(int x=0; x < volSample.xsize();x++)
         if (mask.value(x,y,z) > minValue)
		   {
             sample[i].index = (i+1);
		       sample[i].value = volSample.value(x,y,z);
				 i++;
			}
    sample[i].index = -1;
}

void volume2Sample(svm_model *model, const char *volumeFileName, const char *maskFileName, int sampleSize, float minValue, svm_node * &sample)
{
	volume<float> volSample;
    if (volumeFileName != NULL)
    {
      string volFile = volumeFileName;
      read_volume(volSample, volFile);
    }
   
	volume<float> mask;
    if (maskFileName != NULL)
    {
      string maskFile = maskFileName;
      read_volume(mask, maskFile);
    }
	_volume2Sample(model, volSample, mask, sampleSize, minValue, sample);
}

// function to return the prediction class and svm score of a volume.
float predict(svm_model *model, const char *volumeFileName, const char *maskFileName, int sampleSize, float minValue, float &predictedClass, float &svmScore)
{
	double dPredictedClass, dSvmScore;
	svm_node *sample = NULL;
   volume2Sample(model, volumeFileName, maskFileName, sampleSize, minValue, sample);
	dPredictedClass = svmSamplePredict(model, sample, dSvmScore);
	predictedClass = dPredictedClass;
	svmScore=dSvmScore;
	free(sample);
	return (float) dPredictedClass;
}

// function to return the prediction class and svm score of a volume.
float predict(svm_model *model, const char *volumeFileName, const char *maskFileName, float &predictedClass, float &svmScore)
{
   float minValue=0;
   int sampleSize=getVolumeSize(maskFileName, minValue);
   return predict(model, volumeFileName, maskFileName, sampleSize, minValue, predictedClass, svmScore);
}

// transforms an array in a volume
void array2Volume(const char *maskFile, float minValue, vector <double> &weightVector, volume<float> &weightVolume)
{
	volume<float> mask;
    if (maskFile != NULL)
    {
      string Maskfile = maskFile;
      read_volume(mask, Maskfile);
    }

	weightVolume.reinitialize(mask.xsize(), mask.ysize(), mask.zsize(), 0, true);
   weightVolume.copyproperties(mask);
	int i = 0;
	for(int z=0;z < mask.zsize();z++)
	   for(int y=0;y < mask.ysize();y++)
		  for(int x=0;x < mask.xsize();x++)
			 if (mask.value(x,y,z) > minValue)
			 {
				 weightVolume.value(x,y,z) = (float) weightVector[i];
				 i++;
			 }
			 else weightVolume.value(x,y,z) = (float) 0.0;
}

// the following functions generates a weight volume from a svm model
void _generateWeightVolume(volume<float> &weightVolume, svm_model *model, const char *maskFileName, int vectorSize, float minValue, int normalize)
{
	vector <double> weightVector;
   weightVector.resize(vectorSize);
   
	calculateWeightVector(model, weightVector);
	double norm = 0;
	for (int j=0;j < vectorSize;j++) norm += (weightVector[j] * weightVector[j]);
	if ((norm > 0.0) && (normalize))
	{
		norm = sqrt(norm);
		for (int j=0;j < vectorSize;j++) weightVector[j] /= norm;
	}
	array2Volume(maskFileName, minValue, weightVector, weightVolume);
}

void generateWeightVolume(svm_model *model, const char *maskFileName, int normalize, char *outputFileName)
{
   generateWeightVolume(model, maskFileName, getVolumeSize(maskFileName, 0), 0, normalize, outputFileName);
}

void generateWeightVolume(svm_model *model, const char *maskFileName, int vectorSize, float minValue, int normalize, char *outputFileName)
{
    volume<float>weightVolume;
	_generateWeightVolume(weightVolume, model, maskFileName, vectorSize, minValue, normalize);
    string outname = outputFileName;
    save_volume(weightVolume, outname);
}


