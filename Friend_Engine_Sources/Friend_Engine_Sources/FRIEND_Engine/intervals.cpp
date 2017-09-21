#include <dirent.h>
#include "intervals.h"
#include <string>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <math.h> 
#include <fstream>
#include <list>
#include "defs.h"
#include "fslfuncs.h"
#include "filefuncs.h"
#include <algorithm>


using std::list;

// tests if a condition is the baseline condition
bool DesignObject::isBaselineCondition(const string condition)
{
   return isBaselineCondition(condition.c_str());
}


// tests if a condition is the baseline condition
bool DesignObject::isBaselineCondition(const char *condition)
{
   if (strcmp(baselineCondition, condition) == 0) return 1;
   else return 0;
}

// verifies if the volume index ìdx` is in a baseline block
bool DesignObject::isBaselineCondition(int idx)
{
   return isBaselineCondition(getCondition(idx));
}

// returns the first baseline interval prior to `interval`
int DesignObject::getLastBaseline(int interval)
{
	int baseline = -1;
	if (baselineCondition[0]==0)
	{
		if (interval>0)
		{
         // going back to the first baseline interval before `interval`
		   for (int t=interval-1;t>=0;t--)
		   {
            if (strcmp(baselineCondition, intervals[t].condition)==0)
            {
               baseline=t+1;
               break;
            }
		   }
		}
      // not Found? take the first interval of the run
		if (baseline==-1)
		{
			for (int t=interval;t<intervals.size();t++)
			{
			   if (strcmp(baselineCondition, intervals[t].condition)==0)
			   {
				   baseline=t+1;
				   break;
			   }
			}
		}
	}
	else baseline = interval;
	return baseline;
}

// saves the name of the conditions in a file
void DesignObject::saveConditionNames(char *file)
{
	FILE *f;
	f=fopen(file, "wt+");
	if (f!=NULL)
	{
		for (int t=0;t<conditionNames.size();t++)
      fprintf(f, "%s\n", conditionNames[t].c_str());
		fclose(f);
	}
}

// get condition name of a volume index
const char *DesignObject::getCondition(int idx)
{
   int Idx = returnInterval(idx);
   return (const char*) intervals[Idx].condition;
}

// populates the conditionNames vector with the design file information and sort
int DesignObject::getConditionsList()
{
	conditionNames.clear();
	conditionNamesWithoutBaseline.clear();

   // inserting names that doesn't exists in the list
	for (int i = 0; i< intervals.size(); i++)
	{
      // name doesn't exists. Insert
      if (stringIndex(conditionNames, intervals[i].condition) == -1)
		{
		   string actualCondition = intervals[i].condition;
		   conditionNames.push_back(actualCondition);
		   if (!isBaselineCondition(actualCondition))
			   conditionNamesWithoutBaseline.push_back(actualCondition);
		}
	}
   
   // ordering conditions
   std::sort(conditionNames.begin(), conditionNames.end());
   std::sort(conditionNamesWithoutBaseline.begin(), conditionNamesWithoutBaseline.end());

   return conditionNames.size();
}

// returns the used volume indexes, taking in account offset and baseline conditions
int DesignObject::getVolumeIndices(int offset, int start, int end, vector <int> &indexes, int NoFirstInterval)
{
   indexes.clear();
   
   int initInterval=0;
   if (NoFirstInterval) initInterval=1;
	for (int i=initInterval;i<intervals.size();i++)
	{

		for (int j = (intervals[i].start + offset); j <= intervals[i].end; j++)
		{
			if ((j >= start) && (j <= end) && ((!isBaselineCondition(intervals[i].condition)) || (conditionNamesWithoutBaseline.size()==1))) // if we jut have one condition include baseline
			{
				indexes.push_back(j);
			}
		}
	}
	return indexes.size();
}

// get a condition volume indexes
int DesignObject::getConditionVolumeIndexes(vector <int> &classVector, int classIndex, vector <int> &indexes, vector <int> &classIndexes)
{
   classIndexes.clear();
	for (int t=0;t<indexes.size();t++)
	{
		if (classes[indexes[t]-1]==classIndex)
      classIndexes.push_back(indexes[t]);
	}
	return classIndexes.size();
}

// returns the box car regressor of a condition
int DesignObject::getConditionBoxCar(vector <float> &boxcar, const char*condition, int start, int end)
{
	if (start==0) start = 1;
	if (end==0) end = intervals[intervals.size()-1].end;
	boxcar.clear();
   boxcar.resize(end-start+1);
	for (int i=start;i<=end;i++)
	{
		int interv = returnInterval(i);
		if (interv >= 0)
		{
		   if (strcmp(condition, intervals[interv].condition)==0) boxcar[i-start] = 1;
			else boxcar[i-start]=0;
		}
	}
	return boxcar.size();
}

// save the box car regressor of the conditions in files, one for each condition
void DesignObject::generateConditionsBoxCar(char *outputDir, int start, int end)
{
   int designSize;
   vector <float>boxcar;
   char fileName[BUFF_SIZE];
   
   includeTrailingPathDelimiter(outputDir);
   if (start == 0) start = 1;
   if (end == 0) end = maxIndex();
   for (int i = 0; i < conditionNames.size(); i++)
   {
      getConditionBoxCar(boxcar, conditionNames[i].c_str(), start, end);
      sprintf(fileName, "%s%s%s", outputDir, conditionNames[i].c_str(), ".txt");
      fstream output(fileName, fstream::in | fstream::out | fstream::trunc);
      for (int j=0;j<boxcar.size(); j++) output << boxcar[j] << '\n';
      output.close();
      boxcar.clear();
   }
}

// applies an shift in the array classes. This function will be used for accomodating haemodynamic delays
void DesignObject::shiftClassArray(vector <int> &classesVector, int shift)
{
	for (int t=classesVector.size()-1;t>=0;t--)
	{
		if (t-shift>=0) classesVector[t]=classesVector[t-shift];
		else classesVector[t]=classesVector[0];
	}
}

// returns the condition index in the conditionNames vector of the interval of the given volume index
int DesignObject::getClass(int index)
{
   // get the interval number of the given volume index
   int interval = returnInterval(index);
   // return the condition index of this interval
   return stringIndex(conditionNames, intervals[interval].condition)+1;
}

// populates the class array
int DesignObject::getClassArray(vector<int> &classArray)
{
   classArray.clear();
   classArray.resize(maxIndex(), 0);
	for (int i=0;i<intervals.size();i++)
	{
		int classe = stringIndex(conditionNames, intervals[i].condition)+1;
		for (int j=intervals[i].start;j<=intervals[i].end;j++) classArray[j-1] = classe;
	}
	return classArray.size();
}

// returns the index of name in vector list
int DesignObject::stringIndex(vector <string> &list, char* name)
{
	int r = -1;
	for (int i = 0; i < list.size(); i++)
	{
		if (strcmp(list[i].c_str(), name) == 0)
		{
			r = i;
			break;
		}
	}
	return r;
}

// return the interval of a volue index
int DesignObject::returnInterval(int index)
{
	int ret = -1;
	for (int i=0;i<intervals.size();i++)
	{
		if ((index >= intervals[i].start) && (index <= intervals[i].end))
		{
			ret = i;
			break;
		}
	}
	return ret;
}

// returns the max volume index in the design file
int DesignObject::maxIndex()
{
   return intervals[intervals.size()-1].end;
}

// remove spaces from string
string trim(const std::string &s)
{
	string::const_iterator it = s.begin();
	while (it != s.end() && isspace(*it))
		it++;

	string::const_reverse_iterator rit = s.rbegin();
	while (rit.base() != it && isspace(*rit))
		rit++;

	return string(it, rit.base());
}

// reads a design file
int DesignObject::readDesignFile(char *designFile)
{
	cleanUp();
	ifstream in;
	string row;
	in.open(designFile);
	
	intervals.clear();
   
   // reading a line
	while (getline(in, row))
	{
		row = trim(row);
		if (row.empty()) continue;
	   intervals.resize(intervals.size()+1);
	   {
         // parsing start-end,CONDITION line
         char num[20];
		   int pos = row.find("-");
		   int pos2= row.find(",");
		   row.copy(num, pos);
         
         // start reading
		   num[pos]='\0';
         intervals[intervals.size()-1].start = atoi(num);
		   row.copy(num, pos2-pos-1, pos+1);
         
         // end reading
		   num[pos2-pos-1]='\0';
         intervals[intervals.size()-1].end = atoi(num);
         
         // conditionName reading
		   row.copy(intervals[intervals.size()-1].condition, row.length()-pos2-1, pos2+1);
		   intervals[intervals.size()-1].condition[row.length()-pos2-1] = '\0';
         stripReturns(intervals[intervals.size()-1].condition);
	   }
	}
	in.close();
	getConditionsList();
	return intervals.size();
}

void DesignObject::noInterestContrastsIndexes(char *noInteresetConditions, vector<int> &indexes)
{
	string conditions = trim(string(noInteresetConditions));
	vector<string>names;
	int pos = 0;
	do
	{
		pos = conditions.find(";");
		if (pos > -1)
		{
			string name;
			name = conditions.substr(0, pos);
			names.push_back(name);
			conditions = conditions.substr(pos + 1, conditions.size() - pos - 1);
		}

	} while (pos > -1);
	names.push_back(conditions);

	for (int t = 0; t < contrasts.size(); t++)
	{
		for (int j = 0; j < names.size(); j++)
		{
			if (strstr(contrasts[t].name, names[j].c_str()))
			{
				indexes.push_back(t);
				break;
			}
		}
	}
}

// populates the contrasts and contrast vectors
void DesignObject::generateContrasts(bool conditionContrasts, bool temporalDerivatives)
{
   int contrastInd=0;
   int arrayContrastSize=0, arrayContrastInd=0;
   int contrastSize=0;
   
   // Number of non baseline conditions contrasts 
   contrastSize = conditionNamesWithoutBaseline.size();
   
   // Number of contrasts between non baseline conditions
   if (conditionContrasts)
   {
      // Contrasts with Baseline (T negative condition, J positive condition)
      for (int t=1; t<= conditionNamesWithoutBaseline.size(); t++)
      {
		  for (int j=1; j <= conditionNamesWithoutBaseline.size(); j++)
			  if (t!=j) contrastSize++;
      };
   };
   
   // Defining the size of the contrast vector, depending of the use of temporal derivatives
   if (temporalDerivatives) arrayContrastSize = 2 * conditionNamesWithoutBaseline.size();
   else arrayContrastSize = conditionNamesWithoutBaseline.size();
   
   // Actually filling the contrast information
   // First condition contrast 
   contrasts.resize(contrastSize);
   if (conditionNamesWithoutBaseline.size())
   {
      for (int t = 1; t<= conditionNamesWithoutBaseline.size(); t++)
      {
        contrasts[contrastInd].contrastVector = (float  *) malloc(arrayContrastSize * sizeof(float));
		// zeroing the contrast vector
		for (int j = 0; j < arrayContrastSize; j++) contrasts[contrastInd].contrastVector[j] = 0;

        arrayContrastInd=0;
        for (int j=1; j <= conditionNamesWithoutBaseline.size(); j++)
        {
            if (t==j) contrasts[contrastInd].contrastVector[arrayContrastInd++] = 1;
            else contrasts[contrastInd].contrastVector[arrayContrastInd++] = 0;
            if (temporalDerivatives) contrasts[contrastInd].contrastVector[arrayContrastInd++] = 0;
        };
        sprintf(contrasts[contrastInd].name, "%s", conditionNamesWithoutBaseline[t-1].c_str());
        contrastInd++;
      };
   };
   
   // Now contrast between non baseline conditions
   if (conditionContrasts)
   {
      // Contrasts without Baseline (T negative condition, J positive condition)
      for (int t=1; t <= conditionNamesWithoutBaseline.size(); t++)
      {
        for (int j=1; j <= conditionNamesWithoutBaseline.size(); j++)
        {
            if (t!=j)
            {
                contrasts[contrastInd].contrastVector = (float  *) malloc(arrayContrastSize * sizeof(float));

				// zeroing the contrast vector
				for (int k = 0; k < arrayContrastSize; k++) contrasts[contrastInd].contrastVector[k] = 0;

				arrayContrastInd = 0;
                for (int s=1; s <= conditionNamesWithoutBaseline.size(); s++)
                {
                    if (t==s) contrasts[contrastInd].contrastVector[arrayContrastInd++] = -1;
                    else if (s==j) contrasts[contrastInd].contrastVector[arrayContrastInd++] = 1;
                    else contrasts[contrastInd].contrastVector[arrayContrastInd++] = 0;
                    if (temporalDerivatives) contrasts[contrastInd].contrastVector[arrayContrastInd++] = 0;
                };
                sprintf(contrasts[contrastInd].name, "%s%s%s", conditionNamesWithoutBaseline[j-1].c_str(), "-", conditionNamesWithoutBaseline[t-1].c_str());
                contrastInd++;
            };
        };
      };
   };
};

// generates a simple design matrix. Not used right now
void DesignObject::generateDesignMatrix(char *fileName, int start, int end)
{
   int arraySize;
   vector <float> boxcar;
   vector <vector <float> > convolutions;
   
   convolutions.resize(conditionNames.size());
   
   if (start == 0) start = 1;
   if (end == 0) end = maxIndex();
   arraySize = end-start+1;
   for (int i = 0; i < conditionNames.size(); i++)
   {
      boxcar.clear();
      convolutions[i].resize(arraySize);
      getConditionBoxCar(boxcar, conditionNames[i].c_str(), start, end);
      //convolve(arraySize, 0.72, boxcar, convolutions[i]);
   }
   FILE *f;
   f = fopen(fileName, "wt+");
   char line[255];
   
   for (int i = 0; i < arraySize; i++)
   {
      for (int j = 0; j < conditionNames.size(); j++)
      {
         if (j==0) sprintf(line, "%1.2f", convolutions[j][i]);
         else sprintf(line, "%s\t%1.2f", line, convolutions[j][i]);
      }
      fprintf(f, "%s\n", line);
   }
   fclose(f);
   
   for (int i = 0; i < conditionNames.size(); i++) convolutions[i].clear();
   
   convolutions.clear();
}

// generates the FSF file necessary to run feat_model
void DesignObject::generateFSFFile(char *fileName, int numDynamics, bool motionParams, bool temporalDerivatives)
{
   char customFile[BUFF_SIZE];
   int contrastVectorSize;
   fstream outputFile(fileName, fstream::in | fstream::out | fstream::trunc);

   includeTrailingPathDelimiter(glmDir);
   outputFile << "# Analysis level" << '\n';
   outputFile << "# 1 : First-level analysis" << '\n';
   outputFile << "# 2 : Higher-level analysis" << '\n';
   outputFile << "set fmri(level) 1" << '\n';
   outputFile << "" << '\n';
   outputFile << "# TR(s)" << '\n';
   outputFile << "set fmri(tr) 2" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Total volumes" << '\n';
   outputFile << "set fmri(npts) " << numDynamics << '\n';
   outputFile << "" << '\n';
   outputFile << "# Delete volumes" << '\n';
   outputFile << "set fmri(ndelete) 0" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Critical z for design efficiency calculation" << '\n';
   outputFile << "set fmri(critical_z) 5.3" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Noise level" << '\n';
   outputFile << "set fmri(noise) 0.66" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Noise AR(1)" << '\n';
   outputFile << "set fmri(noisear) 0.34" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Add motion parameters to model" << '\n';
   outputFile << "# 0 : No" << '\n';
   outputFile << "# 1 : Yes" << '\n';
   if (motionParams) outputFile << "set fmri(motionevs) 1" << '\n';
   else outputFile << "set fmri(motionevs) 0" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Number of EVs" << '\n';
   outputFile << "set fmri(evs_orig) " << conditionNamesWithoutBaseline.size() << '\n';
   if (temporalDerivatives) outputFile << "set fmri(evs_real) " << 2*conditionNamesWithoutBaseline.size() << '\n';
   else outputFile << "set fmri(evs_real) " << conditionNamesWithoutBaseline.size()  << '\n';
   outputFile << "set fmri(evs_vox) 0" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Number of contrasts" << '\n';
   outputFile << "set fmri(ncon_orig) " <<  contrasts.size()  << '\n';
   outputFile << "set fmri(ncon_real) " <<  contrasts.size()  << '\n';
   outputFile << "" << '\n';
   outputFile << "# Number of F-tests" << '\n';
   outputFile << "set fmri(nftests_orig) 0" << '\n';
   outputFile << "set fmri(nftests_real) 0" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Highpass temporal filtering" << '\n';
   outputFile << "set fmri(temphp_yn) 1" << '\n';
   outputFile << "" << '\n';
   outputFile << "# Lowpass temporal filtering" << '\n';
   outputFile << "set fmri(templp_yn) 0" << '\n';
   outputFile << "" << '\n';
   outputFile << "# High pass filter cutoff" << '\n';
   outputFile << "set fmri(paradigm_hp) 100" << '\n';
   outputFile << "" << '\n';

   for (int t=1; t<=conditionNamesWithoutBaseline.size(); t++)
   {
      outputFile << "# EV " << t << " title" << '\n';
      outputFile << "set fmri(evtitle" << t << ") \"" << conditionNamesWithoutBaseline[t-1] << "\"" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Basic waveform shape (EV " << t << ")" << '\n';
      outputFile << "# 0 : Square" << '\n';
      outputFile << "# 1 : Sinusoid" << '\n';
      outputFile << "# 2 : Custom (1 entry per volume)" << '\n';
      outputFile << "# 3 : Custom (3 column format)" << '\n';
      outputFile << "# 4 : Interaction" << '\n';
      outputFile << "# 10 : Empty (all zeros)" << '\n';
      outputFile << "set fmri(shape" << t << ") 2" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Convolution (EV " << t << ")" << '\n';
      outputFile << "# 0 : None" << '\n';
      outputFile << "# 1 : Gaussian" << '\n';
      outputFile << "# 2 : Gamma" << '\n';
      outputFile << "# 3 : Double-Gamma HRF" << '\n';
      outputFile << "# 4 : Gamma basis functions" << '\n';
      outputFile << "# 5 : Sine basis functions" << '\n';
      outputFile << "# 6 : FIR basis functions" << '\n';
      outputFile << "set fmri(convolve" << t << ") 3" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Convolve phase (EV " << t << ")" << '\n';
      outputFile << "set fmri(convolve_phase" << t <<") 0" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Apply temporal filtering (EV " << t << ")" << '\n';
      outputFile << "set fmri(tempfilt_yn" << t << ") 1" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Add temporal derivative (EV " << t << ")" << '\n';
      if (temporalDerivatives) outputFile << "set fmri(deriv_yn" << t << ") 1" << '\n';
      else outputFile << "set fmri(deriv_yn" << t << ") 0" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Custom EV file (EV " << t << ")" << '\n';

      sprintf(customFile, "%s%s%s",  glmDir, conditionNamesWithoutBaseline[t-1].c_str(), ".txt");

      outputFile << "set fmri(custom" << t << ") \"" << customFile << "\"" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Gamma sigma (EV " << t << ")" << '\n';
      outputFile << "set fmri(gammasigma" << t << ") 3" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Gamma delay (EV " << t << ")" << '\n';
      outputFile << "set fmri(gammadelay" << t << ") 6" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Orthogonalise EV " << t << " wrt EV 0" << '\n';
      outputFile << "set fmri(ortho" << t << ".0) 0" << '\n';
      for (int j = 1; j <= conditionNamesWithoutBaseline.size(); j++)
      {
         outputFile << "" << '\n';
         outputFile << "# Orthogonalise EV " << t << " wrt EV " << j << '\n';
         outputFile << "set fmri(ortho" << t << "." << j << ") 0" << '\n';
      };
      outputFile << "" << '\n';
   };
   outputFile << "# Contrast & F-tests mode" << '\n';
   outputFile << "# real : control real EVs" << '\n';
   outputFile << "# orig : control original EVs" << '\n';
   outputFile << "set fmri(con_mode_old) orig" << '\n';
   outputFile << "set fmri(con_mode) orig" << '\n';
   outputFile << "" << '\n';

   for (int t = 1; t<= contrasts.size(); t++)
   {
      outputFile << "# Display images for contrast_real " << t << '\n';
      outputFile << "set fmri(conpic_real." << t << ") 1" << '\n';
      outputFile << "" << '\n';
      outputFile << "# Title for contrast_real " << t << '\n';
      outputFile << "set fmri(conname_real." << t << ") \"" << contrasts[t-1].name << "\"" << '\n';
      if (temporalDerivatives) contrastVectorSize = 2 * conditionNamesWithoutBaseline.size();
      else contrastVectorSize = conditionNamesWithoutBaseline.size();
      for (int j = 1; j <= contrastVectorSize; j++)
      {
         outputFile << "" << '\n';
         outputFile << "# Real contrast_real vector " << t << " element " << j << '\n';
         outputFile << "set fmri(con_real" << t << "." << j << ") " << contrasts[t-1].contrastVector[j-1] << '\n';
      };
      outputFile << "" << '\n';
   };
   outputFile.close();
}

// function that cleans up the memory
void DesignObject::cleanUp()
{
   intervals.clear();
   classes.clear();
   conditionNames.clear();
   conditionNamesWithoutBaseline.clear();
   contrasts.clear();
}
