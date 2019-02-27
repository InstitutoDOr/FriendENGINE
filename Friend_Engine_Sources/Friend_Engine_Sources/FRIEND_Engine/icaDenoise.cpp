#include "fMRIProcessing.h"
#include "fMRIUtils.h"

int icaDenoisedSmart(char *iname, char *oname, char *matrix, char *fixIdxs, Matrix& design, int aggressive)
{
	//
	Matrix idxs = read_ascii_matrix(string(fixIdxs));
	Matrix timeSeries = read_ascii_matrix(string(matrix));
	Matrix regressors;
	int valor;


	if (idxs.Nrows())
	{
		znormalization(timeSeries);
		valor = (int)idxs(1, 1);
		regressors = timeSeries.Column(valor);

		for (int ctr = 2; ctr <= (int)idxs.Ncols(); ++ctr){
			valor = (int)idxs(1, ctr);
			regressors |= timeSeries.Column(valor);
		}

		volume4D<float> data;
		read_volume4D(data, string(iname));
		Matrix tmpData = data.matrix();
		Matrix residual;

		//write_ascii_matrix(regressors, "e:\\1.txt");
		regressOut(design, regressors, residual);
		regressors = residual;
		//write_ascii_matrix(regressors, "e:\\2.txt");


		if (aggressive)
		{
			residual = tmpData - (regressors * pinv(regressors) * tmpData);
		}
		else
		{
			Matrix betas = pinv(timeSeries)*tmpData;

			valor = (int)idxs(1, 1);
			Matrix betasSel = betas.Row(valor);

			for (int ctr = 2; ctr <= (int)idxs.Ncols(); ++ctr){
				valor = (int)idxs(1, ctr);
				betasSel &= betas.Row(valor);
			}
			residual = tmpData - regressors * betasSel;
		}

		data.setmatrix(residual);

		save_volume4D(data, string(oname));
	}
	return 0;
}

int icaDenoised(char *iname, char *oname, char *matrix, char *fixIdxs, int aggressive = 0)
{
	//
	Matrix idxs = read_ascii_matrix(string(fixIdxs));
	Matrix timeSeries = read_ascii_matrix(string(matrix));
	Matrix regressors;
	int valor;


	if (idxs.Nrows())
	{
		znormalization(timeSeries);
		valor = (int)idxs(1, 1);
		regressors = timeSeries.Column(valor);

		for (int ctr = 2; ctr <= (int)idxs.Ncols(); ++ctr){
			valor = (int)idxs(1, ctr);
			regressors |= timeSeries.Column(valor);
		}

		volume4D<float> data;
		read_volume4D(data, string(iname));
		Matrix tmpData = data.matrix();
		Matrix residual;


		if (aggressive)
		{
			residual = tmpData - (regressors * pinv(regressors) * tmpData);
		}
		else
		{
			Matrix betas = pinv(timeSeries)*tmpData;

			valor = (int)idxs(1, 1);
			Matrix betasSel = betas.Row(valor);

			for (int ctr = 2; ctr <= (int)idxs.Ncols(); ++ctr){
				valor = (int)idxs(1, ctr);
				betasSel &= betas.Row(valor);
			}
			residual = tmpData - regressors * betasSel;
		}

		data.setmatrix(residual);

		save_volume4D(data, string(oname));
	}
	return 0;
}

int icaResidual(char *iname, char *oname, char *matrix, char *fixIdxs, int aggressive = 0)
{
	//
	Matrix idxs = read_ascii_matrix(string(fixIdxs));
	Matrix timeSeries = read_ascii_matrix(string(matrix));
	Matrix regressors;
	int valor;


	if (idxs.Nrows())
	{
		znormalization(timeSeries);
		valor = (int)idxs(1, 1);
		regressors = timeSeries.Column(valor);

		for (int ctr = 2; ctr <= (int)idxs.Ncols(); ++ctr){
			valor = (int)idxs(1, ctr);
			regressors |= timeSeries.Column(valor);
		}


		volume4D<float> data;
		read_volume4D(data, string(iname));
		Matrix tmpData = data.matrix();
		Matrix residual;

		if (aggressive)
		{
			residual = (regressors * pinv(regressors) * tmpData);
		}
		else
		{
			Matrix betas = pinv(timeSeries)*tmpData;

			valor = (int)idxs(1, 1);
			Matrix betasSel = betas.Row(valor);

			for (int ctr = 2; ctr <= (int)idxs.Ncols(); ++ctr){
				valor = (int)idxs(1, ctr);
				betasSel &= betas.Row(valor);
			}
			residual = regressors * betasSel;
		}

		data.setmatrix(residual);

		save_volume4D(data, string(oname));
	}
	return 0;
}

int icaDenoised(char *dirMelodic, char *fixIdxs, int aggressive = 0)
{
	char iname[2048];
	char oname[2048];
	char matrix[2048];

	sprintf(iname, "%s\\filtered_func_data", dirMelodic);
	sprintf(oname, "%s\\filtered_func_data_cleaned", dirMelodic);
	sprintf(matrix, "%s\\melodic_mix", dirMelodic);

	return icaDenoised(iname, oname, matrix, fixIdxs, aggressive);
}

int icaDenoisedSmart(char *dirMelodic, char* designFile, int aggressive = 0)
{
	char iname[2048];
	char oname[2048];
	char matrix[2048];
	char fixIdxs[2048];
	char onameStd[2048];
	char standard[2048];
	char matrixFile[2048];

	Matrix design, data;
	data = read_ascii_matrix(designFile);
	design = data.Columns(1, 1);

	sprintf(iname, "%s\\filtered_func_data", dirMelodic);
	sprintf(oname, "%s\\filtered_func_data_cleaned", dirMelodic);
	sprintf(matrix, "%s\\filtered_func_data.ica\\melodic_mix", dirMelodic);
	sprintf(fixIdxs, "%s\\TM.txt", dirMelodic);

	int r = icaDenoisedSmart(iname, oname, matrix, fixIdxs, design, aggressive);

	sprintf(onameStd, "%s_std", oname);
	sprintf(standard, "%s\\reg\\standard", dirMelodic);
	sprintf(matrixFile, "%s\\reg\\example_func2standard.mat", dirMelodic);

	//data = read_ascii_matrix(matrixFile);

	//applyTransformation(oname, onameStd, standard, matrixFile);

	return r;
}
