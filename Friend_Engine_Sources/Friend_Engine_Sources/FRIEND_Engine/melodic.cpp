#include "fMRIProcessing.h"

void applySusan(char *entrada, char *saida, char *mascara, float FWHM)
{
	std::stringstream CmdLn;
	volume4D<float> vol;
	volume<float>mask;
	float brain_thresh = 10.0;

	read_volume4D(vol, string(entrada));
	double p2 = vol.percentile(2.0 / 100.0);
	double p98 = vol.percentile(98.0 / 100.0);

	double intensity_threshold = brain_thresh * (p98 - p2) / 100.0;
	fprintf(stderr, "brain thresh = %f p2 = %f p98 = %f int = %f\n", brain_thresh, p2, p98, intensity_threshold);

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << entrada << " -thr " << intensity_threshold << " -Tmin -bin mask -odt char";
	fslmaths((char *)CmdLn.str().c_str());

	read_volume(mask, string("mask"));
	double median_intensity = vol.percentile(50.0 / 100.0, mask);
	fprintf(stderr, "median = %f\n", median_intensity);

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths mask -dilF mask";
	fslmaths((char *)CmdLn.str().c_str());

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << entrada << " -mas mask " << saida;
	fslmaths((char *)CmdLn.str().c_str());


	double smoothsigma = FWHM / 2.355;
	double susan_int = (median_intensity - p2) * 0.75;

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << saida << " -Tmean mean_func";
	fslmaths((char *)CmdLn.str().c_str());

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "susan " << saida << " " << susan_int << "  " << smoothsigma << " 3 1 1 mean_func " << susan_int << " " << saida;
	fprintf(stderr, "%s\n", CmdLn.str().c_str());
	susan((char *)CmdLn.str().c_str());
}

void preparaDataMelodic(char *entrada, char *saida, float FWHM, float highpass, float lowpass, float brain_thresh, float TR, int slicetimerType)
{
	std::stringstream CmdLn;
	volume4D<float> vol;
	volume<float>mask;

	char workFolder[BUFF_SIZE];
	extractFilePath(saida, workFolder);

	char maskFile[BUFF_SIZE];
	char renMaskFile[BUFF_SIZE];

	sprintf(maskFile, "%s\\mean_func.nii.gz", workFolder);
	remove(maskFile);

	sprintf(maskFile, "%s\\mask_mask.nii.gz", workFolder);
	sprintf(renMaskFile, "%s\\mask.nii.gz", workFolder);
	remove(maskFile);
	remove(renMaskFile);

	if (slicetimerType)
	{
		if (slicetimerType == 2) CmdLn << "slicetimer -i " << entrada << " --out=" << saida << " -r " << TR << " --down";
		else if (slicetimerType == 3) CmdLn << "slicetimer -i " << entrada << " --out=" << saida << " -r " << TR << " --odd";
		slicetimer((char *)CmdLn.str().c_str());
		CmdLn.str("");
		entrada = saida;
	}
	CmdLn.precision(10);
	CmdLn << "fslmaths " << entrada << " " << saida << " -odt float";
	fslmaths((char *)CmdLn.str().c_str());

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << saida << " -Tmean " << workFolder << "\\mean_func";
	fslmaths((char *)CmdLn.str().c_str());

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "bet " << workFolder << "\\mean_func " << workFolder << "\\mask -f 0.3 -n -m";
	bet((char *)CmdLn.str().c_str());

	moveFile(maskFile, renMaskFile);

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << saida << " -mas " << workFolder << "\\mask " << saida;
	fslmaths((char *)CmdLn.str().c_str());

	read_volume4D(vol, string(saida));
	double p2 = vol.percentile(2.0 / 100.0);
	double p98 = vol.percentile(98.0 / 100.0);

	double intensity_threshold = brain_thresh * (p98 - p2) / 100.0;
	fprintf(stderr, "brain thresh = %f p2 = %f p98 = %f int = %f\n", brain_thresh, p2, p98, intensity_threshold);

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << saida << " -thr " << intensity_threshold << " -Tmin -bin " << workFolder << "\\mask -odt char";
	fslmaths((char *)CmdLn.str().c_str());

	read_volume(mask, string(renMaskFile));
	double median_intensity = vol.percentile(50.0 / 100.0, mask);
	fprintf(stderr, "median = %f\n", median_intensity);

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << workFolder << "\\mask -dilF " << workFolder << "\\mask";
	fslmaths((char *)CmdLn.str().c_str());

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << saida << " -mas " << workFolder << "\\mask " << saida;
	fslmaths((char *)CmdLn.str().c_str());


	double smoothsigma = FWHM / 2.355;
	double susan_int = (median_intensity - p2) * 0.75;

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << saida << " -Tmean " << workFolder << "\\mean_func";
	fslmaths((char *)CmdLn.str().c_str());

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "susan " << saida << " " << susan_int << "  " << smoothsigma << " 3 1 1 " << workFolder << "\\mean_func " << susan_int << " " << saida;
	fprintf(stderr, "%s\n", CmdLn.str().c_str());
	susan((char *)CmdLn.str().c_str());


	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << saida << " -mas " << workFolder << "\\mask " << saida;
	fslmaths((char *)CmdLn.str().c_str());

	double normmean = 10000.0;
	double scaling = normmean / median_intensity;

	CmdLn.str("");
	CmdLn.precision(10);
	CmdLn << "fslmaths " << saida << " -mul " << scaling << " " << saida;
	fprintf(stderr, "%s\n", CmdLn.str().c_str());
	fslmaths((char *)CmdLn.str().c_str());

	if ((highpass > 0) || (lowpass > 0))
	{
		float hp_sigma_vol = -1;
		if (highpass > 0)
		{
			float hp_sigma_sec = highpass / 2.0;
			hp_sigma_vol = hp_sigma_sec / TR;
		}

		double lp_sigma_vol = -1;
		if (lowpass > 0)
		{
			double lp_sigma_sec = 2.8;
			lp_sigma_vol = lp_sigma_sec / TR;
		}

		CmdLn.str("");
		CmdLn.precision(10);
		CmdLn << "fslmaths " << saida << " -bptf " << hp_sigma_vol << " " << lp_sigma_vol << " " << saida;
		fprintf(stderr, "%s\n", CmdLn.str().c_str());
		fslmaths((char *)CmdLn.str().c_str());
	}
}

void runMelodic(char *dirSuj, char *dirRun)
{
	char dirMelodic[BUFF_SIZE];
	char inName[BUFF_SIZE], outName[BUFF_SIZE];

	sprintf(dirMelodic, "%s\\Melodic", dirSuj);
	sprintf(inName, "%s\\sl%s.nii", dirSuj, dirRun);

	mkdir(dirMelodic);

	sprintf(outName, "%s\\sl%s", dirMelodic, dirRun);

	preparaDataMelodic(inName, outName);

	char cmd[1024];
	sprintf(cmd, "melodic -i %s -o %s.ica -v --nobet --bgthreshold=3 --tr=2.0 --report --guireport=%s\\report.html -d 0 --mmthresh=0.5 --bgimage=%s\\mean_func --vn --Ostats", outName, outName, dirMelodic, dirMelodic);
	melodic(cmd);
}

void runMelodic4D(char *inName)
{
	char dirMelodic[BUFF_SIZE];
	char outName[BUFF_SIZE];

	changeFileExt(inName, ".ica", dirMelodic);

	mkdir(dirMelodic);

	sprintf(outName, "%s\\filtered_func_data.nii", dirMelodic);

	preparaDataMelodic(inName, outName, 6.0, 100.0, 0);

	char cmd[1024];
	sprintf(cmd, "melodic -i %s -o %s.ica -v --nobet --bgthreshold=3 --tr=2.0 --report --guireport=%s\\report.html -d 0 --mmthresh=0.5 --bgimage=%s\\mean_func --vn --Ostats", outName, outName, dirMelodic, dirMelodic);
	melodic(cmd);
}
