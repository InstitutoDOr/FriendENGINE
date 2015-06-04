//
//  masks.h
//  
//
//  Created by Radiologia Rededor on 5/21/13.
//
//

#ifndef ____masks__
#define ____masks__

#include <iostream>
#include <vector>
#include <map>
#include "newimage/newimageall.h"

using namespace NEWIMAGE;
using namespace std;

// struct for roi processing
typedef struct
{
   int x, y, z;
   float value;
   int roiValue;
   
} roiPoint;

// structs for vector class sorting
typedef struct
{
   bool operator()(roiPoint const &a, roiPoint const &b) const { return a.value > b.value; }
} greaterRoiPoint;

typedef struct
{
   bool operator()(roiPoint const &a, roiPoint const &b) const { return a.value < b.value; }
} lesserRoiPoint;

// znormalize a vetor of double
void znormalise(vector<double>&values);

// calculates the mean of a set of rois
class RoiMeanCalculation
{
private:
   map<int, int> mapping;
   vector<float> means;
   vector<int> counts;
   volume<int> reference;
   
public:
	// loads the volume that defines a set of rois based on the voxels intensities. For that works the voxel datatype must be int
   void loadReference(char *referenceFileName);
   // initializes the object variables
   void reinitialize();
   // calculate the rois means of a given volume
   void calculateMeans(char *volumeFileName);
   // calculates the mean for each roi
   void calculateMeans(volume<float> &actualvolume);
   // returns the voxels values of reference volume 
   int voxelValues(int x, int y, int z);
   // returns the index of a intensity roi  
   int roiIndex(int roiValue);
   // returns the number of voxels in a roi  
   int roiSize(int roiIndex);
   // returns the mean value of a voxel  
   float roiMean(int roiIndex);
   // returns the number of rois
   int roiCount();
   // returns the number of rois
   void getRoiValues(vector<int> &values);
};

// extracts the best voxels of each roi
class RegionExtraction
{
   RoiMeanCalculation reference;
public:
   // returns in output the best voxels of `region`
   void regionBestVoxels(RoiMeanCalculation &reference, volume<float>&values, volume<float>&output, int region, int regionSize, float percentage);
   // extract the `percentage` best voxels of each roi, based on a volume that defines the rois and volumes of T value. How these volumes are used is explained in code
   void regionsExtraction(char *refVol, char *valueVol4D, char *valueVol, char *outputVol, map<int,int> &regionContrastMap, float percentage);

   // reaads a region extraction map file
   void readMappings(char *filename, std::map<int, int> &mappings);
   
};

// calculates mean, standard deviation and variance incrementally
class IncrementalStats
{
public:
    float mean, variance, deviation;
    float minValue, maxValue;
    int numPts, numPtsVariance;
    float window;
   
    // update the statistical variables with another value, modifying the variance if required (default is modify)
    void addValue(float value, bool modifyvariance = 1);
   
    // checks if a value is inside a deviation window from mean
    bool valueInsideDeviationWindow(float value, float devMult=1);
   
    // value inside a slack window from mean
    bool valueInsideWindow(float value);
   
    // calculates the zscore based on these incremental variables
    double zValue(double value);
   
    // initializing variables
    void initialize();
};

// calculates a weighted mean a list of values incrementally increased
class WeightedMean
{
public:
   vector<float> vectorData;
   vector<float> meanCoefs;
   int width;
   IncrementalStats incStats;
   float window;
   float mean;
   int usedPoints;
   
   // initializing variables
   void initialize();
   
   // sets how many values are used in mean calculation
   void setMeanWidth(int value);
   
   // sets an initial list of values size
   void setVectorSize(int size);
   
   // value inside a slack window from mean
   bool valueInsideWindow(float value);
   
   // calculates the mean
   float meanExtract(int index);
   
   // initialize mean coeficients
   void initializeCoefs(int type = 1);
   
   // updates the state of variables with another value
   void addValue(float value, int calculateMean = 1);

   // saves the curve and the mean curve in a file, separated by ;
   void saveCurves(char *file);
   
#ifndef __GNUC__
   WeightedMean() {  width = 0; };
#else //The GCC way
   WeightedMean() { width = 0; };
#endif
   virtual ~WeightedMean() { };
};

typedef struct
{
    int map, roiIntensity;
}
regionRoi;


typedef struct
{
    int regionA, regionB;
}
regionPair;


// calculates the correlation between rois
class RegionCorrelation {
public:
    int numRegions;
    int numPairs;
    int runSize;
    int finals;
    int correlationType;
    int windowSize;
    int actualSize;
    bool cachedMovingWindow;
    
    vector<RoiMeanCalculation> calculatorMaps;
    vector< vector<float> > regionVector;
    vector<regionPair> pairs;
    vector< vector<regionRoi> > regions;

    void _calculateMeans();
    void _calculateCorreations(int idx);
    double correlation(vector<float>A, vector<float>B, int nFinals);
   
    
public:
    // vector of region means
    vector<float> means;
   // vector of correlations, one for each pair of regions
    vector<float> correlations;

    // sets the size of the list of roi maps
    void setCalculatorMapsSize(int size);
    // loads a roi map in one place in the map list
    void loadVOIS(int map, char *filename);
    // this is the most basic and only used right now. Just one map in the list
    void loadRegionMap(char *filename);
    // adds a new region to the list of regions. Remember that one region can be composed of more than one roi, defined by its roi's intensity in the roi map
    void addRegion(vector<int> &rois, int map = 1);
    // adds a pair of regions in correlation pair list. For now we have only one pair
    void addPair(int regionA, int regionB);
    // sets the size of region list, for further assigning of values
    void setRegionNumber(int numRegion);
    // returns the size of roi map list
    int  getMapNumber();
    // sets the size of correlation pairs list, for further assigning of region values
    void setPairNumber(int pairs);
    void loadRegions(char *filename, char *referenceDirectory = NULL);
    // zeroes the region mean vector
    void zeroVectors();
    // saves the actual configuration. Note : for now it does not save the map list
    void saveRegions(char *filename);
    // calculates the composed roi region mean of a volume
    void calculateMeans(char *filename);
    // restarting the calculation. Actual size holds the last presented volume number. See next function for details
    void resetCalculation();
    void addMeanDataPoint(int idx, float value);
    // sets the correlation window size
    void setWindowSize(int value);
    // sets the region vector size
    void setRunSize(int value);
    // updates the correlation calculation with a new volume
    void calculateCorrelations(char *volume, int idx=0);
    // updates the correlation calculation with a new volume
    void calculateCorrelations(volume<float> &volcorr, int idx=0);

    
    RegionCorrelation()
    {
        finals = 0;
        numRegions = 0;
        runSize = 0;
        correlationType = 1;
    }
};

#endif /* defined(____masks__) */
