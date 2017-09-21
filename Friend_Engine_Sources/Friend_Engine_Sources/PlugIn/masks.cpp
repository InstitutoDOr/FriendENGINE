//
//  masks.cpp
//  
//
//  Created by Radiologia Rededor on 5/21/13.
//
//

#include "statistics.h"
#include "ap.h"
#include <cmath>
#include "masks.h"
#include <fstream>
#include "parser.h"
#include "filefuncs.h"
#include "defs.h"

using namespace alglib;

// loads the volume that defines the initial values
void EMADetrend::loadReference(char *referenceFileName, float alphaValue)
{
	volume<float>vol;
	load_volume(vol, string(referenceFileName));
	loadReference(vol, alphaValue);
}

// loads the volume that defines the initial values
void EMADetrend::loadReference(volume<float>&volume, float alphaValue)
{
	emaFiltered.reinitialize(volume.xsize(), volume.ysize(), volume.zsize());
	emaFiltered.copyproperties(volume);

	int numVoxels = emaFiltered.xsize() * emaFiltered.ysize() * emaFiltered.zsize();
	float *meanDataPtr = (float *)emaFiltered.fbegin();
	float *volDataPtr = (float *)volume.fbegin();
	for (int t = 0; t < numVoxels; t++) meanDataPtr[t] = volDataPtr[t];

	alpha = alphaValue;
	numPoints = 1;
}

// update the baseline calculation
void EMADetrend::addVolume(char *volumeName)
{
	volume<float>vol;
	load_volume(vol, string(volumeName));
	addVolume(vol);
}

// update the baseline calculation
void EMADetrend::addVolume(volume<float> &vol)
{
	numPoints++;
	int numVoxels = emaFiltered.xsize() * emaFiltered.ysize() * emaFiltered.zsize();
	float *meanDataPtr = (float *)emaFiltered.fbegin();
	float *volDataPtr = (float *)vol.fbegin();

	float emaFilt, value;
	for (int t = 0; t < numVoxels; t++)
	{
		value = volDataPtr[t];
		emaFilt = meanDataPtr[t];
		emaFilt = emaFilt * alpha + (1 - alpha) * value;
		meanDataPtr[t] = emaFilt;
		volDataPtr[t] = value - emaFilt;
	}
}

// loads a value that defines the initial value
void IncrementalBaseline::loadReference(float value)
{
	numVoxels = 1;
	meanArray = (baselineStruct *)malloc(numVoxels * sizeof(baselineStruct));

	meanArray[0].mean = value;
	meanArray[0].numPoints = 1;
}

// loads the volume that defines the initial values
void IncrementalBaseline::loadReference(char *referenceFileName, int zeroValues)
{
	volume<float>vol;
	load_volume(vol, string(referenceFileName));
	loadReference(vol, zeroValues);
}

// loads the volume that defines the initial values
void IncrementalBaseline::loadReference(volume<float>&volume, int zeroValues)
{
	numVoxels = volume.xsize() * volume.ysize() * volume.zsize();
	meanArray = (baselineStruct *) malloc(numVoxels * sizeof(baselineStruct));

	float *volDataPtr = (float *) volume.fbegin();
	for (int t = 0; t < numVoxels; t++)
	{
		if (zeroValues) meanArray[t].mean = 0;
		else meanArray[t].mean = volDataPtr[t];
		meanArray[t].numPoints = 1;
	}
}

// update the baseline calculation
void IncrementalBaseline::addVolume(char *volumeName, int subtract)
{
	volume<float>vol;
	load_volume(vol, string(volumeName));
	addVolume(vol, subtract);
}

// update the baseline calculation
void IncrementalBaseline::addVolume(volume<float> &vol, int subtract)
{
	float *volDataPtr = (float *)vol.fbegin();

	float mean, oldMean, value;
	for (int t = 0; t < numVoxels; t++)
	{
		value = volDataPtr[t];
		mean = meanArray[t].mean;
		if (mean > value)
		{
			meanArray[t].numPoints++;
			oldMean = mean;
			mean = mean + (value - oldMean) / meanArray[t].numPoints;
			meanArray[t].mean = mean;
		}
		if (subtract)
			volDataPtr[t] = value - mean;
	}
}

// update the baseline calculation
void IncrementalBaseline::addVolume(float &value, int subtract)
{
	float mean, oldMean;
	mean = meanArray[0].mean;
	if (mean > value)
	{
		meanArray[0].numPoints++;
		oldMean = mean;
		mean = mean + (value - oldMean) / meanArray[0].numPoints;
		meanArray[0].mean = mean;
	}
	if (subtract)
		value = value - mean;
}

// znormalize a vetor of double
void znormalise(vector<double>&values)
{
	IncrementalStats stats;
	stats.initialize();
	for (int i = 0; i < values.size(); i++) stats.addValue(values[i]);
	for (int i = 0; i < values.size(); i++) values[i] = stats.zValue(values[i]);
}

// creates a vector of weights with points of a sigmoid fuction
void sigmoid(int dataPoints, vector<float> &points, float slope = 1, float start = -6, float end = 6)
{
    float incr, value;
    incr = (end-start)/(dataPoints-1);
    points.resize(dataPoints);
    value = start;
    for (int  t=0;t < dataPoints;t++)
    {
        points[t] = 1/(1 + exp(-slope*value));
        value += incr;
    }
}

// returns the index of a voxel intensity in an vector of roi intensities
bool insideRegion(float value, vector<int>regions)
{
   int intValue = (int) value;
   bool result=false;
   for (int t=0; t < regions.size();t++)
      if (regions[t]==intValue)
      {
         result=true;
         break;
      }
   return result;
}

// loads the volume that defines a set of rois based on the voxels intensities. For that works the voxel datatype must be int
void RoiMeanCalculation::loadReference(char *referenceFileName)
{
   reinitialize();
   if (!fileExists(referenceFileName))
   {
	   fprintf(stderr, "Reference volume %s not found.\n", referenceFileName);
	   return;
   }
   read_volume(reference, string(referenceFileName));
   
   // counts the number of rois
   int idx=0;
   for(int z=0;z<reference.zsize();z++)
      for(int y=0;y<reference.ysize();y++)
         for(int x=0;x<reference.xsize();x++)
         {
            int value=reference.value(x,y,z);
            if (value != 0)
            {
               std::pair<std::map<int,int>::iterator,bool> ret = mapping.insert( std::pair<int, int>(value,idx++) );
               // alredy exists. Decrementing the count
               if (ret.second==false) idx--;
            }
         }
   
   // initiate the mean and count vectors
   means.resize(idx);
   counts.resize(idx);
   for (int t=0; t<idx; t++)
   {
      means[t]=0;
      counts[t]=0;
   }
   
   // Calculate the sie of each roi
   for(int z=0;z<reference.zsize();z++)
      for(int y=0;y<reference.ysize();y++)
         for(int x=0;x<reference.xsize();x++)
         {
            int value=reference.value(x,y,z);
            if (value != 0)
            {
               // here the mapping object returns the associated index paired with the intensity value
               int index = mapping[value];
               counts[index]++;
            }
         }
}

// initializes the object variables
void RoiMeanCalculation::reinitialize()
{
   mapping.clear();
   means.clear();
   counts.clear();
}

// calculate the rois means of a given volume
void RoiMeanCalculation::calculateMeans(char *volumeFileName)
{
   volume<float> actualVolume;
   if (means.size())
   {
	   read_volume(actualVolume, string(volumeFileName));
	   calculateMeans(actualVolume);
   }
}

// calculates the mean for each roi
void RoiMeanCalculation::calculateMeans(volume<float> &actualvolume)
{
	if (means.size())
	{
		for (int i = 0; i < means.size(); i++) means[i] = 0;

		for (int z = 0; z < reference.zsize(); z++)
			for (int y = 0; y < reference.ysize(); y++)
				for (int x = 0; x < reference.xsize(); x++)
				{
					float value = reference.value(x, y, z);
					if (value != 0)
					{
						int idx = mapping[value];
						if (idx >= -1)
							means[idx] += actualvolume.value(x, y, z);
					}
				}

		for (int i = 0; i < means.size(); i++)
		{
			if (counts[i]) means[i] = means[i] / (float)counts[i];
		}
	}
}

// returns the voxels values of reference volume 
int RoiMeanCalculation::voxelValues(int x, int y, int z)
{
	return reference.value(x, y, z);
}

// returns the index of a intensity roi  
int RoiMeanCalculation::roiIndex(int roiValue)
{
	return mapping[roiValue];
}

// returns the number of voxels in a roi  
int RoiMeanCalculation::roiSize(int roiIndex)
{
	if (roiIndex < counts.size())
	{
		return counts[roiIndex];
	}
	else return 0;
}

// returns the number of rois
int RoiMeanCalculation::roiCount()
{
	return means.size();
}

// returns the mean value of a voxel  
float RoiMeanCalculation::roiMean(int index)
{
	if (index < means.size())
	{
		return means[index];
	}
	else return 0;
}

// returns the number of rois
void RoiMeanCalculation::getRoiValues(vector<int> &values)
{
	values.clear();
	for (std::map<int, int>::iterator it = mapping.begin(); it != mapping.end(); ++it)
		values.push_back(it->first);
}

// returns the roi value given its index
int RoiMeanCalculation::getRoiValue(int index)
{
	vector<int> values;
	getRoiValues(values);

	if (index < values.size())
	{
		return values[index];
	}
	else return 0;
}

// reads a region extraction map file
void RegionExtraction::readMappings(char *filename, std::map<int, int> &mappings)
{
   char buffer[500];
   int count;
   char **numbers;
   
   fstream input(filename, fstream::in | fstream::out);
   if (input.is_open())
   {
      while (input.good())
      {
         input.getline(buffer, 100);
         // parsing the line in the for <roi intensity> = <contrast number>
         parser(buffer, count, numbers, '=');
         if (count>1)
            mappings[strtol(numbers[0], NULL, 10)]=strtol(numbers[1], NULL, 10);
         freeparser(count, numbers);
      }
   }
}

// returns in `output` the best voxels of `region`
void RegionExtraction::regionBestVoxels(RoiMeanCalculation &reference, volume<float>&values, volume<float>&output, int region, int regionSize, float percentage)
{
   vector<roiPoint> roi;
   roi.resize(regionSize);
   greaterRoiPoint greaterFirst;
   int index=0;
   
   // Filter the voxels of the region `region`
   for(int z=0;z<values.zsize();z++)
      for(int y=0;y<values.ysize();y++)
         for(int x=0;x<values.xsize();x++)
         {
            // get the voxel intensity in reference
            int voxelRegion = (int) reference.voxelValues(x,y,z);
            
            // if is the chosen region, records the voxel values (T values or other voxel value)
            if (region == voxelRegion)
            {
               roi[index].value = values.value(x,y,z);
               roi[index].roiValue = region;
               
               roi[index].x=x;
               roi[index].y=y;
               roi[index].z=z;
               index++;
            }
            else values.value(x,y,z)=(float)0.0;
         }
   
   // sorts the vector of values in descending order
   std::sort(roi.begin(), roi.end(), greaterFirst);
   
   // calculates the cut Index. Remember the voxels descending order
   int cutIndex = (int) (roi.size() * percentage + 0.5);
   
   // recording the result in `output`
   for (int j=0; j<cutIndex;j++)
   {
      output.value(roi[j].x, roi[j].y, roi[j].z) = roi[j].roiValue;
   }
}

// extract the `percentage` best voxels of each roi, based on a volume that defines the rois and volumes of T value. How these volumes are used is explained in code
void RegionExtraction::regionsExtraction(char *refVol, char *valueVol4D, char *valueVol, char *outputVol, map<int,int>&regionContrastMap, float percentage)
{
   // loads the reference volume demarking the rois
   reference.loadReference(refVol);
   vector< vector<roiPoint> >rois;
   vector<int> indexes;
   
   volume<float>values;
   volume4D<float>values4D;
   volume<float>output;
   
   // loads the summary volume
   read_volume(values, string(valueVol));
   
   // loads the 4D volume of contrasts
   read_volume4D(values4D, string(valueVol4D));
   
   // initializating the volume output variable
   output.reinitialize(values.xsize(), values.ysize(), values.zsize());
   output.copyproperties(values);
   
   // this is important. Zeroing the volume
   output=0;

   // calculate each region separately. regionContrastMap defines the map of a roi intensity and a contrast index.
   for (std::map<int,int>::iterator it = regionContrastMap.begin(); it!=regionContrastMap.end(); ++it)
   {
      // getting region <-> contrast duo
      int region = it->first;
      int contrast = it->second;
      
      // gets the region index, to get the size of the region
      int idxregion = reference.roiIndex(region);
      
      // if contrast > 0, get the values from the volume contrast
      if (contrast) regionBestVoxels(reference, values4D[contrast-1], output, region, reference.roiSize(idxregion), percentage);
      else
         // get the values from the summary volume
		 regionBestVoxels(reference, values, output, region, reference.roiSize(idxregion), percentage);
   }
   // saving the result
   save_volume(output, string(outputVol));
}

// extract the `percentage` best voxels of a roi, based on a volume that defines the roi and T value volume
void RegionExtraction::regionExtraction(char *refVol, char *valueVol, char *outputVol, float percentage)
{
	// loads the reference volume demarking the roi
	reference.loadReference(refVol);

	volume<float>values;
	volume<float>output;

	// loads the summary volume
	read_volume(values, string(valueVol));

	// initializating the volume output variable
	output.reinitialize(values.xsize(), values.ysize(), values.zsize());
	output.copyproperties(values);

	// this is important. Zeroing the volume
	output = 0;

	// calculate each region separately. regionContrastMap defines the map of a roi intensity and a contrast index.
	regionBestVoxels(reference, values, output, reference.getRoiValue(0), reference.roiSize(0), percentage);
	// saving the result
	save_volume(output, string(outputVol));
}

// update the statistical variables with another value, modifying the variance if required (default is modify)
void IncrementalStats::addValue(float value, bool modifyVariance)
{
   float oldMean;
   numPts++;
   
   if (numPts == 1)
   {
      minValue = value;
      maxValue = value;
   }
   else
   {
      minValue = ((minValue > value) ? value : minValue);
      maxValue = ((maxValue < value) ? value : maxValue);
   };
   
   oldMean = mean;
   mean = mean + (value-oldMean)/numPts;
   
   if (modifyVariance)
   {
      numPtsVariance++;
      variance += (value-oldMean)*(value-mean);
      if (numPtsVariance < 2) deviation = 0;
      else deviation = sqrt(variance/(numPtsVariance-1));
   }
}

// checks if a value is inside a deviation window from mean
bool IncrementalStats::valueInsideDeviationWindow(float value, float devMult)
{
    return (value >= mean - deviation * devMult) && (value <= mean + deviation * devMult);
}

// value inside a window from mean
bool IncrementalStats::valueInsideWindow(float value)
{
    return (value >= mean - window) && (value <= mean + window);
}

// calculates the zscore of a value based on these incremental variables
double IncrementalStats::zValue(double value)
{
    if (deviation == 0) return 0;
    else return (value - mean) / deviation;
}

// initializing variables
void IncrementalStats::initialize()
{
    mean = 0;
    variance = 0;
    deviation = 0;
    numPts = 0;
    numPtsVariance = 0;
    minValue = 0;
    maxValue = 0;
    window = 0;
};

// initializing variables
void WeightedMean::initialize()
{
    usedPoints = 0;
    mean = 0;
    window = 0;
};

// sets how many values are used in mean calculation
void WeightedMean::setMeanWidth(int value)
{
    meanCoefs.resize(value);
    width = value;
    sigmoid(value, meanCoefs);
};

// sets an initial list of values size
void WeightedMean::setVectorSize(int size)
{
    usedPoints = 0;
    vectorData.resize(size);
};

// value inside a slack window from mean
bool WeightedMean::valueInsideWindow(float value)
{
  return (value >= mean - window) && (value <= mean + window);
};

// calculates the mean. Index starts with one
float WeightedMean::meanExtract(int index)
{
    int startIndex;
    float coefsSum;
   
    // initializing variables
    startIndex = index - width -1;
    startIndex = (startIndex > 0) ? startIndex : 0;
   
    mean = 0;
    coefsSum = 0;
    incStats.initialize();
    // this calculates the weighted mean from the last value added backwards the first in the width
    for (int i = index-1; i >= startIndex; i--)
    {
        // backward weighting
        mean = mean + vectorData[i] * meanCoefs[(width-1)+i-(index-1)];
       
        // this for another calculation not used anymore
        incStats.addValue(vectorData[i]);
       
        // constructing the denominator
        coefsSum = coefsSum + meanCoefs[(width-1)+i-(index-1)];
    }

	// returning the result
	if (coefsSum != 0) mean = mean / coefsSum;
	else mean = 0;
    return mean;
};

// initialize mean coefficients
void WeightedMean::initializeCoefs(int type)
{
   if (type == 1) sigmoid(meanCoefs.size(), meanCoefs);
   else
   for(int t=0; t < meanCoefs.size(); t++) meanCoefs[t] = 1;
};

void WeightedMean::calculateLimits(int lastValues)
{
	int first = usedPoints - lastValues;
	if (first < 0) first = 0;
	for (int i = first; i < usedPoints; i++)
	{
		if (i == first)
		{
			minValue = vectorData[i];
			maxValue = vectorData[i];
		}
		else
		{
			if (minValue > vectorData[i]) minValue = vectorData[i];
			if (maxValue < vectorData[i]) maxValue = vectorData[i];
		}
	};
}

// converts value to a percentage inside the min-max range. 
// Use the slack variable if you want to avoid ceiling effect
float WeightedMean::scaleValue(float value, float slack)
{
	float tempMax = maxValue, tempMin = minValue;
	if (slack)
	{
		tempMax *= (1 + slack);
		tempMin *= (1 - slack);
	};

	float range = tempMax - tempMin;
	if (range == 0) return 0;
	else return (value - tempMin) / (range);
}

// converts value to a percentage inside the min-max range. 
// Use the slack variable if you want to avoid ceiling effect
// This function calls the calculate limits to resolve things in one call
float WeightedMean::scaleValue(float value, float slack, int lastValues)
{
	calculateLimits(lastValues);
	return scaleValue(value, slack);
}

// updates the state of variables with another value
void WeightedMean::addValue(float value, int calculateMean)
{
   // if no more room, resize the vector. Its important to preset the size correctly to avoid a lot of resizings
   if (vectorData.size() < usedPoints+1) vectorData.resize(usedPoints+1);
   
   vectorData[usedPoints] = value;
   usedPoints++;
   
   // calculating the mean
   if (calculateMean) meanExtract(usedPoints);
};

// saves the curve and the mean curve in a file, separated by ;
void WeightedMean::saveCurves(char *file)
{
	fstream output(file, fstream::in | fstream::out | fstream::trunc);
	for (int i = 0; i < usedPoints; i++)
	{
		meanExtract(i+1);
		output << vectorData[i] << ";" << mean << "\n";
	}
	output.close();
}

// initialize the variables
void RegionCorrelation::initialize()
{
	numRegions=0;
	numPairs=0;
	runSize=0;
	finals=0;
	correlationType=1;
	windowSize=10;
	actualSize=0;
	cachedMovingWindow=true;

	calculatorMaps.resize(0);
	regionVector.resize(0);
	regions.resize(0);

	pairs.resize(0);
	// vector of region means
	means.resize(0);
	correlations.resize(0);

}

// sets the size of the list of roi maps
void RegionCorrelation::setCalculatorMapsSize(int size)
{
   calculatorMaps.resize(size);
};

// loads a roi map in one place in the roi map list
void RegionCorrelation::loadVOIS(int map, char *filename)
{
   calculatorMaps[map-1].loadReference(filename);
};

// this is the most basic and only used right now. Just one map in the list
void RegionCorrelation::loadRegionMap(char *filename)
{
   setCalculatorMapsSize(1);
   loadVOIS(1, filename);
};

// adds a new region to the list of regions. Remember that one region can be composed of more than one roi, defined by its roi's intensity in the roi map
void RegionCorrelation::addRegion(vector<int> &rois, int map)
{
   if (rois.size())
   {
      // increasing the proper variables
      numRegions++;
      regions.resize(numRegions);
      means.resize(numRegions);

      regionVector.resize(numRegions);
      regionVector[numRegions-1].resize(runSize);
      
      // making room for the intensity roi list of that region
      regions[numRegions-1].resize(rois.size());
      
      // setting the list of roi intensities of the region
      for (int t=0; t < rois.size(); t++)
      {
         // by default, map = 1 (just one roi map)
         regions[numRegions-1][t].map = map;
         regions[numRegions-1][t].roiIntensity = rois[t];
      }
   }
}

// adds a pair of regions in correlation pair list. For now we have only one pair
void RegionCorrelation::addPair(int regionA, int regionB)
{
   numPairs++;
   
   pairs.resize(numPairs);
   pairs[numPairs-1].regionA = regionA;
   pairs[numPairs-1].regionB = regionB;

   correlations.resize(numPairs);
   correlations[numPairs-1] = 0;
}

// sets the size of region list, for further assigning of values
void RegionCorrelation::setRegionNumber(int numRegion)
{
   numRegions = numRegion;
   regions.resize(numRegions);
   means.resize(numRegions);

   regionVector.resize(numRegions);
   for (int t = 0; t<numRegions; t++)
   {
      regionVector[t].resize(runSize);
   }
}

// returns the size of roi map list
int RegionCorrelation::getMapNumber()
{
   return calculatorMaps.size();
}

// sets the size of correlation pairs list, for further assigning of region values
void RegionCorrelation::setPairNumber(int pairsNumber)
{
   numPairs = pairsNumber;
   pairs.resize(numPairs);
   
   correlations.resize(numPairs);
   for (int t=0; t<numPairs;t++) correlations[t]=0;
}

// reads a roi region definition. It can be of the form region{-map}, optionally indicating from wich map we read the roi information. For now, the map indication is not used.
void getRoiInfo(char* roidef, int &roiIntensity, int &map)
{
   map = 1;
   char *divider = strchr(roidef, '-');
   if (divider)
   {
      // there is a map indication
      divider[0]='\0';
      roiIntensity = atoi(roidef);
      divider++;
      map = atoi(divider);
   }
   else roiIntensity = atoi(roidef);
}

// reads a region file definition of regions, pirs of correlations nd possibly roi map files
void RegionCorrelation::loadRegions(char *filename, char *referenceDirectory)
{
   fstream regionFile(filename, std::fstream::in);
   char bufferLine[512];
   
   regionFile.getline(bufferLine, sizeof(bufferLine)-1);
   
   // reading regions. First line is the number of regions
   int regionNumber = atoi(bufferLine);
   if (regionNumber)
   {
      setRegionNumber(regionNumber);
      for (int t=0; t<regionNumber;t++)
      {
         int roisNumber;
         char **rois;

         regionFile.getline(bufferLine, sizeof(bufferLine)-1);
         parser(bufferLine, roisNumber, rois, ',');
         regions[t].resize(roisNumber);
         for (int i=0;i<roisNumber;i++)
         {
            int roiIntensity, map;
            
            getRoiInfo(rois[i], roiIntensity, map);
            regions[t][i].roiIntensity = roiIntensity;
            regions[t][i].map = map;
         }
         freeparser(roisNumber, rois);
      }
   }
   
   // now reading correlation pairs
   regionFile.getline(bufferLine, sizeof(bufferLine)-1);

   // first line is the number of pairs
   int pairNumber = atoi(bufferLine);
   setPairNumber(pairNumber);
   
   for (int t=0; t<pairNumber; t++)
   {
      int size;
      char **charPairs;
      
      regionFile.getline(bufferLine, sizeof(bufferLine)-1);
      parser(bufferLine, size, charPairs, ',');
      if (size ==2)
      {
         pairs[t].regionA = atoi(charPairs[0]);
         pairs[t].regionB = atoi(charPairs[1]);
      }
      freeparser(size, charPairs);
   }
   
   // reading the maps list, if any.
   regionFile.getline(bufferLine, sizeof(bufferLine)-1);
   if (bufferLine[0]==0)
   {
      // First line is the number of maps
      int numMaps = atoi(bufferLine);
      for (int t=0;t<numMaps; t++)
      {
         char mapFile[BUFF_SIZE];
         // reading the map file name
         regionFile.getline(bufferLine, sizeof(bufferLine)-1);
         // if not found , add the referenceDir in the name
         if (!fileExists(bufferLine)) sprintf(mapFile, "%s%c%s", referenceDirectory, PATHSEPCHAR, bufferLine);
         else strcpy(mapFile, bufferLine);
         // load the map in list
         if (fileExists(mapFile)) loadVOIS(t+1, mapFile);
      }
   }
   
}

// the most simply way to load information in this object
void RegionCorrelation::loadRoiMask(char *roiMask)
{
	RoiMeanCalculation mask;
	vector<int>rois, values;

	if (fileExists(roiMask))
	{ 
		mask.loadReference(roiMask);

		mask.getRoiValues(values);
		// at least two different rois should exist
		if (values.size() > 1)
		{
			rois.resize(1);

			// first roi value 
			rois[0] = values[0];
			addRegion(rois);

			// second roi value 
			rois[0] = values[1];
			addRegion(rois);

			addPair(1, 2);
			loadRegionMap(roiMask);
		}
	}
}

// zeroes the region mean vector
void RegionCorrelation::zeroVectors()
{
   for (int t=0; t<numRegions; t++)
      for (int j=0; j<runSize; j++)
          regionVector[t][j] = 0;
}

// saves the actual configuration. Note : for now it does not save the map list
void RegionCorrelation::saveRegions(char *filename)
{
   fstream regionFile(filename, std::fstream::out);
    // Regions
    if (numRegions)
    {
        regionFile << numRegions << "\n";
        for (int t=0;t<numRegions; t++)
        {
            for (int j=0; j<regions[t].size(); j++)
            {
                regionFile << regions[t][j].roiIntensity << "-" << regions[t][j].map;
                if (j<regions.size()-1) regionFile << ",";
            }
            regionFile << "\n";
        }
        regionFile << "\n";
    }
    
    // Pairs
    if (numPairs)
    {
        regionFile << numPairs << "\n";
        for (int t=0; t< numPairs; t++)
        {
            regionFile << pairs[t].regionA << "," << pairs[t].regionB << "\n";
        }
    }
}

// calculates the composed roi region mean
void RegionCorrelation::_calculateMeans()
{
    for (int t=0;t<numRegions;t++)
    {
        double mean=0;
        for (int j=0; j<regions[t].size();j++)
        {
            // getting the index of a roi in the region
            int idx = calculatorMaps[regions[t][j].map-1].roiIndex(regions[t][j].roiIntensity);
            // adding the roi mean to the calculation of region mean
            if (idx > -1)
            {
                mean += calculatorMaps[regions[t][j].map-1].roiMean(idx);
            }
        }
       
        // finishing the calculation
        if (regions[t].size()) mean /= regions[t].size();
        else mean = 0;
       
        // and recording the result
        means[t] = mean;
    }
    
}

// calculates the composed roi region mean of a volume
void RegionCorrelation::calculateMeans(char *filename)
{
    // here we are only calculating the roi means
    for(int t=0;t<calculatorMaps.size();t++)
        calculatorMaps[t].calculateMeans(filename);
   
    // here we are actually calculating the region mean
    _calculateMeans();
}

// sets the region vector size
void RegionCorrelation::setRunSize(int value)
{
   runSize=value;
}

// sets the correlation window size
void RegionCorrelation::setWindowSize(int value)
{
    resetCalculation();
    cachedMovingWindow = true;
    windowSize = value;
    if (runSize==0)
    {
       runSize = value;
       for (int t=0; t<numRegions; t++) regionVector[t].resize(value);
    }
    finals = value;
}

// updates the correlation calculation with a new volume
void RegionCorrelation::calculateCorrelations(volume<float>&volcorr, int idx)
{
   for (int t=0;t<calculatorMaps.size();t++) calculatorMaps[t].calculateMeans(volcorr);
   // _calculateMeans() is called in _calculateCorreations
   _calculateCorreations(idx);
}

// updates the correlation calculation with a new volume
void RegionCorrelation::calculateCorrelations(char *volume, int idx)
{
    for (int t=0;t<calculatorMaps.size();t++) calculatorMaps[t].calculateMeans(volume);
    // _calculateMeans() is called in _calculateCorreations
    _calculateCorreations(idx);
}

// calculates the window correlation between to vectors
double RegionCorrelation::correlation(vector<float>A, vector<float>B, int nFinals)
{
   int iniIdx = 0, endIdx = A.size();
   int size;
   if (nFinals)
   {
      iniIdx = endIdx-nFinals; // here we are not adding 1 to account for zero based index
      if (iniIdx < 0) iniIdx = 0;
   };
   // not always size = nFinals
   size = endIdx-iniIdx;
   
   real_1d_array a, b;
   
   // transferring the values to another data structure for correlation calculation (alglib)
   a.setlength(size);
   for (int t=iniIdx;t<endIdx; t++) a[t] = A[t];
   
   b.setlength(size);
   for (int t=iniIdx;t<endIdx; t++) b[t] = B[t];
   
   // here we work with to types of correlation. The defautl is pearson
   if (correlationType == 1) return pearsoncorr2(a, b, size);
   else return spearmancorr2(a, b, size);
}

// restarting the calculation in cached mode. Actual size holds the last presented index in cache. See next function for details
void RegionCorrelation::resetCalculation()
{
   actualSize = 0;
}

// adds a mean data in cache
void RegionCorrelation::addMeanDataPoint(int idx, float value)
{
   if (actualSize > windowSize)
   { // cache full, get rid of the unused stuff
      for (int t=0;t<(windowSize-1);t++)
      {
         regionVector[idx][t] = regionVector[idx][t+1];
      }
   }
   // adding the new value
   regionVector[idx][((actualSize > windowSize) ? windowSize : actualSize)-1] = value;
}

// actually calculates the correlations
void RegionCorrelation::_calculateCorreations(int idx)
{
   if (calculatorMaps.size())
   {
       // calculate composed region means
       _calculateMeans();
      
       // updating the next cache index accordly
       if (cachedMovingWindow)
       {
          actualSize++;
          actualSize = (actualSize > (windowSize+1)) ? windowSize+1 : actualSize;
       }
      
       // adding the result in the vectors
       for (int t=0;t<numRegions;t++)
       {
           if (cachedMovingWindow) addMeanDataPoint(t, means[t]);
           else regionVector[t][idx] = means[t];
       }
       
       // calculate correlations
       for (int t=0; t<numPairs; t++)
       {
           double value;
           // calculating the correlation between a pair of region means
           if (cachedMovingWindow)
           {
               if (windowSize > actualSize) value = 0;
			   else value = correlation(regionVector[pairs[t].regionA - 1], regionVector[pairs[t].regionB - 1], windowSize);
           }
           else value = correlation(regionVector[pairs[t].regionA-1], regionVector[pairs[t].regionB-1], finals);

           // recording the result
           correlations[t] = value;
       }
   }
}
