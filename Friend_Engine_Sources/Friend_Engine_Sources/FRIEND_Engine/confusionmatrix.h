// Functions responsible of confusion matrix calculations

#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

class ConfusionMatrix
{
   public:
   vector<int> matrix;
   vector<string> classNames;
   int matrixRank;
   
   // initializes the confusion matrix
   ConfusionMatrix();
   
   // sets the number of classes of the experiment
   void ConfusionMatrix::setRank(int rank);

   // sets a class name
   void setClassName(int classNumber, const string &className);
   
   // resets the accuracy information of the matrix
   void resetMatrix();
   
   // copy the information of a confusion matrix
   void copyMatrix(const ConfusionMatrix& source);
   
   // making the = operator works
   const ConfusionMatrix& operator=(const ConfusionMatrix& source);

   // add the result information of a sample test
   void reportResult(int target, int acquired);
   
   // return the kappa value of the samples informed to the confusion matrix
   // kappa = (Hits percentage - Chance percentage) / (100 - Chance percentage)
   float kappa();
   
   // calculates the chance probability
   float chance();
   
   // returns the accucary of the samples informed to the confusion matrix
   float hits();
   
   // calculates the number of examples correct classified of a class informed to the confusion matrix
   int examples();

   // calculate the sensibility of a class informed to the confusion matrix
   float sensibilityClassPrediction(int classNumber);
   
   // calculates the number of examples correct classified informed to the confusion matrix
   float accuracyClassPrediction(int classNumber);
   
   // save a text file report of the confusion matrix
   void saveMatrixReport(const char *file);
};
