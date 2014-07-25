#include "confusionmatrix.h"

// initializes the confusion matrix
ConfusionMatrix::ConfusionMatrix(int rank)
{
   matrixRank = rank;
	matrix.resize(matrixRank * matrixRank);
   classNames.resize(matrixRank);
}

// sets a class name. Note : the first class number is 1
void ConfusionMatrix::setClassName(int classNumber, const string &className)
{
   classNames[classNumber-1] = className;
}

// resets the accuracy information of the matrix
void ConfusionMatrix::resetMatrix()
{
	int maximum = matrixRank*matrixRank;
	for (int i=0; i<maximum;i++) matrix[i] = 0;
}

// copy the information of a confusion matrix
void ConfusionMatrix::copyMatrix(const ConfusionMatrix& source)
{
   matrixRank=source.matrixRank;
	int maximum = matrixRank*matrixRank;
   matrix.resize(maximum);
	for (int i=0; i<maximum;i++) matrix[i] = source.matrix[i];
}

// making the = operator works
const ConfusionMatrix& ConfusionMatrix::operator=(const ConfusionMatrix& source)
{
   this->copyMatrix(source);
   return *this;
}

// add the result information of a sample test. Note the target and acquired class numbers starts in 1
void ConfusionMatrix::reportResult(int target, int acquired)
{
	matrix[(target-1) * matrixRank + (acquired-1)]++;
}

// return the kappa value of the samples informed to the confusion matrix
float ConfusionMatrix::kappa() // (Hits percentage - Chance percentage) / (100 - Chance percentage)
{
   return (hits() - chance()) / (100.0-chance());
}

// calculates the chance probability
float ConfusionMatrix::chance()
{
	return (float) (100.0 / (float) matrixRank);
}

// returns the accucary of the samples informed to the confusion matrix
float ConfusionMatrix::hits()
{
   int total=0, acc=0;

   for (int target=0; target < matrixRank; target++)
	   for (int acquired=0; acquired < matrixRank; acquired++)
	   {
		   if (target==acquired) acc = acc + matrix[target*matrixRank + acquired];
		   total = total + matrix[target*matrixRank + acquired];
	   };
   return ((float)acc * 100) / (float) total;
}

// calculates the number of examples informed to the confusion matrix
int ConfusionMatrix::examples()
{
	int maximum = matrixRank*matrixRank;
	int total=0;
	for (int i=0; i<maximum;i++) total+=matrix[i];
	return total;
}

// calculates the number of examples correct classified of a class informed to the confusion matrix
float ConfusionMatrix::accuracyClassPrediction(int classNumber)
{
   int acc=0, total=0;
   for (int i=0;i<matrixRank;i++)
   {
	   if (i == classNumber) acc = matrix[i*matrixRank+classNumber-1];
	   total += matrix[i*matrixRank+classNumber-1];
   }
   return (float)(acc * 100) / (float)total;
}

// calculate the sensibility of a class informed to the confusion matrix
float ConfusionMatrix::sensibilityClassPrediction(int classNumber)
{
   int acc=0, total=0;
   for (int i=0;i<matrixRank;i++)
   {
      if (i == classNumber) acc = matrix[(classNumber-1)*matrixRank+i];
	  total += matrix[(classNumber-1)*matrixRank+i];
   }
   return (float)(acc * 100) / (float)total;
}

// save a text file report of the confusion matrix
void ConfusionMatrix::saveMatrixReport(const char *file)
{
	FILE *f=fopen(file, "wt+");
	char buffer[255];
	if (f!=NULL)
	{
		fprintf(f, "-------------------------------------------------------------------------------------------------------------------------------\n");
        fprintf(f, "|%20s|", "Names");
		for (int i=0;i<matrixRank;i++) fprintf(f, "%20s|", classNames[i].c_str());
        fprintf(f, "%20s|\n", "Accuracy");
		fprintf(f, "-------------------------------------------------------------------------------------------------------------------------------\n");

		for (int i=0;i<matrixRank;i++) 
		{
			fprintf(f, "|%20s|", classNames[i].c_str());
			int acc=0, total=0;
			for (int j=0;j<matrixRank;j++)
			{
				fprintf(f, "%20d|", matrix[i*matrixRank+j]);
				if (i==j) acc = matrix[i*matrixRank+j];
				total += matrix[i*matrixRank+j];

			}
			if (total == 0) sprintf(buffer, "%.2f", 0.0f);
			else sprintf(buffer, "%.2f", (float)acc*100/(float)total);
			strcat(buffer, "%");
			fprintf(f, "%20s|\n", buffer);
		}

		fprintf(f, "-------------------------------------------------------------------------------------------------------------------------------\n");
		fprintf(f, "|%20s|", "Specificity");

		for (int j=0;j<matrixRank;j++)
		{
			int acc=0, total=0;
			for (int i=0;i<matrixRank;i++)
			{
				if (i==j) acc = matrix[i*matrixRank+j];
				total += matrix[i*matrixRank+j];
			}
			if (total == 0) sprintf(buffer, "%.2f", 0.0f);
			else sprintf(buffer, "%.2f", (float)acc*100/(float)total);
			strcat(buffer, "%");
			fprintf(f, "%20s|", buffer);
		}
        fprintf(f, "%20s|\n", "");
		fprintf(f, "-------------------------------------------------------------------------------------------------------------------------------\n");
		fprintf(f, "\n");

		sprintf(buffer, "%.2f", hits());
		strcat(buffer, "%");

		fprintf(f, "kappa Value   : %.2f\n", kappa());
		fprintf(f, "Total Accuracy : %s\n", buffer);
		fprintf(f, "Total Examples : %d\n", examples());
		fclose(f);
	}
}

