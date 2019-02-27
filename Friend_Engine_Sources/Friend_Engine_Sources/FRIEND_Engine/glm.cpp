#include "fMRIProcessing.h"
#include "libprob.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/t2z.h"
#include "miscmaths/f2z.h"


void basicGLM(const Matrix& data, const Matrix& design, const Matrix& contrasts, Matrix &t, int DOFadjust = 0)
{
	Matrix beta;
	Matrix residu;
	Matrix sigsq;
	Matrix varcb;
	Matrix cbeta;
	Matrix f_fmf, pf_fmf;
	int dof;
	Matrix z;
	Matrix p;

	beta = zeros(design.Ncols(), 1);
	residu = zeros(1); sigsq = -1.0*ones(1); varcb = -1.0*ones(1);
	t = zeros(1); z = zeros(1); p = -1.0*ones(1);
	dof = (int)-1; cbeta = -1.0*ones(1);

	if (data.Nrows() == design.Nrows()){
		Matrix dat = data;
		Matrix tmp = design.t()*design;
		Matrix pinvdes = tmp.i()*design.t();

		beta = pinvdes * dat;
		residu = dat - design*beta;

		dof = design.Nrows() - design.Ncols() - 1;
		sigsq = sum(SP(residu, residu)) / dof;

		float fact = float(dof) / design.Ncols();
		f_fmf = SP(sum(SP(design*beta, design*beta)), pow(sum(SP(residu, residu)), -1)) * fact;

		pf_fmf = f_fmf.Row(1);
		for (int ctr1 = 1; ctr1 <= f_fmf.Ncols(); ctr1++)
			pf_fmf(1, ctr1) = 1.0 - MISCMATHS::fdtr(design.Ncols(),
			int(design.Nrows() - 1 - design.Ncols()), f_fmf.Column(ctr1).AsScalar());

		if (contrasts.Storage()>0 && contrasts.Ncols() == beta.Nrows()){
			cbeta = contrasts*beta;
			Matrix tmp = contrasts*pinvdes*pinvdes.t()*contrasts.t();
			varcb = diag(tmp)*sigsq;
			t = SP(cbeta, pow(varcb, -0.5));
			z = t; p = t;
			for (int ctr1 = 1; ctr1 <= t.Ncols(); ctr1++){
				ColumnVector tmp = t.Column(ctr1);
				T2z::ComputeZStats(varcb.Column(ctr1), cbeta.Column(ctr1), dof, tmp);
				z.Column(ctr1) << tmp;
				T2z::ComputePs(varcb.Column(ctr1), cbeta.Column(ctr1), dof, tmp);
				p.Column(ctr1) << exp(tmp);
			}
		}
	}

}


void subjectGLM(char *dirSuj, char *dirRun, int tipo = 1)
{
	char dirGLM[BUFF_SIZE];
	char inName[BUFF_SIZE], outName[BUFF_SIZE], designFile[BUFF_SIZE], contrastFile[BUFF_SIZE];
	volume4D<float>data;

	if (tipo == 1)
	{
		sprintf(inName, "%s\\sl%s.nii", dirSuj, dirRun);
		sprintf(outName, "%s\\sl%s_t.nii", dirSuj, dirRun);
	}
	else if (tipo == 2)
	{
		sprintf(inName, "%s\\ema%s.nii", dirSuj, dirRun);
		sprintf(outName, "%s\\ema%s_t.nii", dirSuj, dirRun);
	}
	else if (tipo == 3)
	{
		sprintf(inName, "%s\\%s.nii", dirSuj, dirRun);
		sprintf(outName, "%s\\%s_t.nii", dirSuj, dirRun);
	}

	read_volume4D(data, string(inName));


	sprintf(dirGLM, "%s\\GLM", dirSuj);


	Matrix tmpData, design, contrasts, t;
	tmpData = data.matrix();
	tmpData = MISCMATHS::remmean(tmpData, 1);

	sprintf(designFile, "%s\\%s_MP.mat", dirGLM, dirRun);
	sprintf(contrastFile, "%s\\%s_MP.con", dirGLM, dirRun);

	design = read_ascii_matrix(string(designFile));
	contrasts = read_ascii_matrix(string(contrastFile));

	basicGLM(tmpData, design, contrasts, t);

	volume4D<float> volumeT(data.xsize(), data.ysize(), data.zsize(), t.Nrows());
	volumeT.setmatrix(t);
	volumeT.copyproperties(data[0]);
	save_volume4D(volumeT, string(outName));
}
