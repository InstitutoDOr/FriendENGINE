#include "fMRIProcessing.h"
#include "fMRIUtils.h"

void regressOut(Matrix &regressors, Matrix&data, Matrix &residual)
{
	int cols = regressors.Ncols();
	int rows = regressors.Nrows();
	Matrix auxReg(rows, cols + 1);

	for (int i = 1; i <= rows; i++)
		for (int j = 1; j <= cols + 1; j++)
			if (j <= cols) auxReg(i, j) = regressors(i, j);
			else auxReg(i, j) = 1;

	Matrix tmp = auxReg.t()*auxReg;
	Matrix pinvdes = tmp.i()*auxReg.t();
	Matrix betaAux = pinvdes * data;

	if (1)
	{
		Matrix beta = betaAux.Rows(1, cols);
		residual = data - regressors*beta;
	}
	else residual = data - auxReg*betaAux;
}

void regressOutDesign(char *inFile, char *outFile, int initialColumn)
{
	Matrix residual;
	Matrix data = read_ascii_matrix(string(inFile));

	Matrix design = data.Columns(1, initialColumn);
	Matrix regressors = data.Columns(initialColumn + 1, data.Ncols());
	regressOut(design, regressors, residual);
	write_ascii_matrix(residual, string(outFile));
}

void regressOutDesign(char *inFile, char *designFile, char *outFile)
{
	Matrix residual;
	Matrix data = read_ascii_matrix(string(inFile));
	Matrix design = read_ascii_matrix(string(designFile));

	regressOut(design, data, residual);
	write_ascii_matrix(residual, string(outFile));
}

void estimateActivation4D(volume4D<float> &data, int slidingWindow)
{
	// iterates for each volume
	int ini = 1;
	int end = data.tsize();
	volume4D<float>activation(data.xsize(), data.ysize(), data.zsize(), end - ini + 1);
	activation.copyproperties(data);

	for (int i = ini; i <= end; i++)
	{
		int beginning = Max(1, i - slidingWindow + 1);
		activation[i - ini] = 0;

		// sliding window mean
		for (int j = beginning; j <= i; j++)
			activation[i - ini] += data[j - 1];
		activation[i - ini] /= (float)(i - beginning + 1);
	}
	data = activation;
}

void regressOutData(volume4D<float> &data, char *regressorFile)
{
	Matrix tmpData = data.matrix();
	Matrix regressors = read_ascii_matrix(string(regressorFile));
	Matrix residual;

	regressOut(regressors, tmpData, residual);
	data.setmatrix(residual);
}

void regressOutFromData(char *iname, char *regressorFile, char *oname)
{
	volume4D<float> data;

	read_volume4D(data, string(iname));
	regressOutData(data, regressorFile);
	save_volume4D(data, string(oname));
}

void regressOutFromDataTrends(char *iname, char *regressorFile, char *oname, int NOrdem)
{
	volume4D<float> data;

	read_volume4D(data, string(iname));
	Matrix tmpData = data.matrix();
	Matrix regressors = read_ascii_matrix(string(regressorFile));
	Matrix extraColumns(regressors.Nrows(), 1);
	Matrix residual;
	for (int n = 1; n <= NOrdem; n++)
	{
		for (int l = 1; l <= regressors.Nrows(); l++)
		{
			if (n == 1) extraColumns(l, 1) = l;
			else extraColumns(l, 1) = pow(l, n);
		}
		regressors |= extraColumns;
	}
	znormalization(regressors);
	regressOut(regressors, tmpData, residual);
	data.setmatrix(residual);
	save_volume4D(data, string(oname));
}