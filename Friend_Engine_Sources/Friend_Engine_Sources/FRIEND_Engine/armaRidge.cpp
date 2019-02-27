#include "armaRidge.h"
#include <cmath>


template<class T> inline T sqr(const T& x) { return x*x; };

void Matrix2Mat(Matrix &a, mat &b)
{
	b.set_size(a.Nrows(), a.Ncols());
	for (int i = 0; i < b.n_rows; i++)
	{
		for (int j = 0; j < b.n_cols; j++) b(i, j) = a(i+1, j+1);
	}
}

void Mat2Matrix(mat &a, Matrix &b)
{
	b.ReSize(a.n_rows, a.n_cols);
	for (int i = 0; i < a.n_rows; i++)
	{
		for (int j = 0; j < a.n_cols; j++) b(i+1, j+1) = a(i, j);
	}
}

void znormalization_arma(mat& matrix, const int dim)
{
	mat matmean = mean(matrix, dim);
	mat matstdev = stddev(matrix, 0, dim);

	if (dim == 0){
		for (int mr = 0; mr < matrix.n_rows; mr++)
			matrix.row(mr) -= matmean.row(0);

		for (int mc = 0; mc < matrix.n_cols; mc++)
			matrix.col(mc) /= matstdev(0, mc);
	}
	else{
		for (int mc = 0; mc < matrix.n_cols; mc++)
			matrix.col(mc) -= matmean.col(0);

		for (int mr = 0; mr < matrix.n_rows; mr++)
			matrix.row(mr) /= matstdev(mr, 0);
	}
}

void ridgeArma(mat &regressors, mat &data, doubleVect &alphas, mat &betas, int normalpha, double singcutoff)
{
	vec S;
	mat U, V;


	svd_econ(U, S, V, regressors, "both");
	float tol;
	tol = max(abs(min(S)), max(S)) * max(regressors.n_rows, regressors.n_cols) * singcutoff;

	double S0;
	int SOSetted = 0, k = 0;
	for (int n = 0; n < S.n_rows; n++) {
		if (fabs(S(n)) > tol)
		{
			if (!SOSetted)
			{
				S0 = S(n);
				SOSetted = 1;
			}
			k++;
		}
		else S(n) = 0.0; // reduce the number of columns because too close to singular
	}

	mat UR = U.t() * data;
	doubleVect nalphas;
	for (int t = 0; t < alphas.size(); t++)
	{
		if (normalpha) nalphas.push_back(alphas[t] * S0);
		else nalphas.push_back(alphas[t]);
	}

	betas.set_size(regressors.n_cols, nalphas.size());

	mat D;
	D = diagmat(S);
	for (int t = 0; t < nalphas.size(); t++)
	{
		for (int n = 0; n < S.n_rows; n++) {
			if (fabs(S(n)) > 0)
			{
				D(n, n) = S(n) / (sqr(S(n)) + nalphas[t]);
			}
			else D(n, n) = 0;
		}
		mat w = V * D * UR;
		for (int r = 0; r < betas.n_rows; r++) betas(r, t) = w(r, 0);
	}
}

void ridgeArma2(mat &regressors, mat &data, doubleVect &alphas, mat &betas, mat &corr, int normalpha, double singcutoff)
{
	vec S;
	mat U, V;

	svd_econ(U, S, V, regressors, "both");
	float tol;
	tol = max(abs(min(S)), max(S)) * max(regressors.n_rows, regressors.n_cols) * singcutoff;

	double S0;
	int SOSetted = 0;
	for (int n = 0; n < S.n_rows; n++) {
		if (fabs(S(n)) > tol)
		{
			if (!SOSetted)
			{
				S0 = S(n);
				SOSetted = 1;
			}
		}
		else S(n) = 0.0; // reduce the number of columns because too close to singular
	}

	doubleVect nalphas;
	for (int t = 0; t < alphas.size(); t++)
	{
		if (normalpha) nalphas.push_back(alphas[t] * S0);
		else nalphas.push_back(alphas[t]);
	}

	mat UR = U.t() * data;
	betas.set_size(regressors.n_cols, nalphas.size());

	mat D, SD;
	SD = diagmat(S);
	D = SD;
	corr.set_size(nalphas.size(), 1);
	for (int t = 0; t < nalphas.size(); t++)
	{
		for (int n = 0; n < S.n_rows; n++) {
			if (fabs(S(n)) > 0)
			{
				D(n, n) = S(n) / (sqr(S(n)) + nalphas[t]);
			}
			else D(n, n) = 0;
		}
		mat DUR = D * UR;
		mat w = V * DUR;
		mat pred = U * SD * DUR;

		znormalization_arma(pred);
		for (int i = 0; i < pred.n_cols; i++)
		{
			double sum = 0;
			for (int j = 0; j < pred.n_rows; j++)
				sum += pred(j, i) * data(j, i);

			corr(t, 0) = sum / (pred.n_rows - 1);
		}
		for (int r = 0; r < betas.n_rows; r++) betas(r, t) = w(r, 0);
	}
}

void ridge_corr(mat &regressorsTrain, mat &regressorsTest, mat &dataTrain, mat &dataTest, doubleVect &alphas, mat &correlations, int normalpha, double singcutoff)
{
	vec S;
	mat U, V;


	svd_econ(U, S, V, regressorsTrain, "both");
	float tol;
	tol = max(abs(min(S)), max(S)) * max(regressorsTrain.n_rows, regressorsTrain.n_cols) * singcutoff;


	double S0;
	int SOSetted = 0;
	for (int n = 0; n < S.n_rows; n++) {
		if (fabs(S(n)) > tol)
		{
			if (!SOSetted)
			{
				S0 = S(n);
				SOSetted = 1;
			}
		}
		else S(n) = 0.0; // reduce the number of columns because too close to singular
	}

	doubleVect nalphas;
	for (int t = 0; t < alphas.size(); t++)
	{
		if (normalpha) nalphas.push_back(alphas[t] * S0);
		else nalphas.push_back(alphas[t]);
	}

	mat UR = U.t() * dataTrain;
	mat PV = regressorsTest * V;

	mat zdata = dataTest;
	znormalization_arma(zdata);

	mat D;
	D = diagmat(S);

	correlations.set_size(nalphas.size(), dataTrain.n_cols);
	for (int t = 0; t < nalphas.size(); t++)
	{
		for (int n = 0; n < S.n_rows; n++) {
			if (fabs(S(n)) > 0)
			{
				D(n, n) = sqr(S(n)) / (sqr(S(n)) + nalphas[t]);
			}
			else D(n, n) = 0;
		}
		mat pred = PV * D * UR;

		znormalization_arma(pred);

		for (int i = 0; i < pred.n_cols; i++)
		{
			double sum = 0;
			for (int j = 0; j < pred.n_rows; j++)
				sum += pred(j, i) * zdata(j, i);

			correlations(t, i) = sum / (pred.n_rows - 1);
		}
	}
}

void calculateBetas(mat &regressors, mat &data, doubleVect &alphas, mat &betas, int normalpha, double singcutoff)
{
	vec S;
	mat U, V;

	svd_econ(U, S, V, regressors, "both");
	float tol;
	tol = max(abs(min(S)), max(S)) * max(regressors.n_rows, regressors.n_cols) * singcutoff;

	double S0;
	int SOSetted = 0;
	for (int n = 0; n < S.n_rows; n++) {
		if (fabs(S(n)) > tol)
		{
			if (!SOSetted)
			{
				S0 = S(n);
				SOSetted = 1;
			}
		}
		else S(n) = 0.0; // reduce the number of columns because too close to singular
	}

	doubleVect nalphas;
	for (int t = 0; t < alphas.size(); t++)
	{
		if (normalpha) nalphas.push_back(alphas[t] * S0);
		else nalphas.push_back(alphas[t]);
	}

	betas.set_size(regressors.n_cols, data.n_cols);

	mat D, SD;
	SD = diagmat(S);
	D = SD;
	for (int t = 0; t < nalphas.size(); t++)
	{
		for (int n = 0; n < S.n_rows; n++) {
			if (fabs(S(n)) > 0)
			{
				D(n, n) = S(n) / (sqr(S(n)) + nalphas[t]);
			}
			else D(n, n) = 0;
		}
		mat w = V * D * U.t() * data.col(t);
		for (int r = 0; r < betas.n_rows; r++) betas(r, t) = w(r, 0);
	}
}

void crossValidation(mat & regressors, mat &data, doubleVect &alphas, doubleVect &bestAlphas, int nfolds, unsigned int seed)
{
	vector<int>kfolds;
	for (int i = 0; i < regressors.n_rows; i++) kfolds.push_back((i % nfolds) + 1);

	bestAlphas.resize(data.n_cols);
	if (seed) srand(seed);
	random_shuffle(kfolds.begin(), kfolds.end());
	mat correlations, meanCorrelations;
	meanCorrelations.set_size(alphas.size(), data.n_cols);
	meanCorrelations.fill(0);

	for (int t = 1; t <= nfolds; t++)
	{
		mat regressorsTrain, regressorsTest, dataTrain, dataTest;
		int numTest = 0;
		for (int i = 0; i < kfolds.size(); i++)
			if (kfolds[i] == t) numTest++;

		regressorsTrain.set_size(kfolds.size()-numTest, regressors.n_cols);
		regressorsTest.set_size(numTest, regressors.n_cols);
		dataTrain.set_size(kfolds.size() - numTest, data.n_cols);
		dataTest.set_size(numTest, data.n_cols);

		int idxTrain=0, idxTest=0;
		for (int i = 0; i < kfolds.size(); i++)
		{
			if (kfolds[i] == t)
			{
				regressorsTest.row(idxTest) = regressors.row(i);
				dataTest.row(idxTest++) = data.row(i);
			}
			else
			{
				regressorsTrain.row(idxTrain) = regressors.row(i);
				dataTrain.row(idxTrain++) = data.row(i);
			}
		}
		fprintf(stderr, "Fold = %d Training = %d Test = %d\n", t, kfolds.size() - numTest, numTest);

		ridge_corr(regressorsTrain, regressorsTest, dataTrain, dataTest, alphas, correlations);
		meanCorrelations += correlations;
	}

	meanCorrelations /= nfolds;
	for (int i = 0; i < meanCorrelations.n_cols; i++)
	{
		double bestAlpha=0, bestCorrelation=-1.0;
		for (int j = 0; j < meanCorrelations.n_rows; j++)
		{
			if (bestCorrelation < meanCorrelations(j, i))
			{
				bestCorrelation = meanCorrelations(j, i);
				bestAlpha = alphas[j];
			}
		}
		bestAlphas[i] = bestAlpha;
	}
}