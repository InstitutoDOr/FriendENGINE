#include "fMRIProcessing.h"
#include "intervals.h"
#include "stringUtils.h"

void filtroGaussiano(char *entrada, char *saida, float FWHM)
{
	char CmdLn[1024];
	sprintf(CmdLn, "fslmaths %s -kernel gauss %f -fmean %s", entrada, FWHM / 2.3548, saida);
	fslmaths(CmdLn);
}

void processData(volume4D<float> &data, float FWHM, char *regressorFile, DesignObject &interval, char *dirSuj, char *dirRun)
{
	volume<float> vol;
	int withBaseline = 0;
	char outputFile[BUFF_SIZE];

	// apagando antigos
	sprintf(outputFile, "%s\\ema%s.nii", dirSuj, dirRun);
	remove(outputFile);
	sprintf(outputFile, "%s\\sl%s.nii", dirSuj, dirRun);
	remove(outputFile);
	sprintf(outputFile, "%s\\%s.nii", dirSuj, dirRun);
	remove(outputFile);

	// Regressores
	fprintf(stderr, "Regressando dados\n");
	regressOutData(data, regressorFile);

	// Sliding Window
	fprintf(stderr, "Sliding Window\n");
	estimateActivation4D(data, 3);

	// Save slidingWindow
	sprintf(outputFile, "%s\\sl%s.nii", dirSuj, dirRun);
	save_volume4D(data, string(outputFile));
	return;

	// Detrending
	fprintf(stderr, "Ema detrending\n");

	// EMA Detrending
	EMADetrend emaObj;
	for (int t = 0; t < data.tsize(); t++)
	{
		if (t == 0)
		{
			emaObj.loadReference(data[t]);
			data[t] = 0;
		}
		else emaObj.addVolume(data[t]);
	}

	// Save EMA
	sprintf(outputFile, "%s\\ema%s.nii", dirSuj, dirRun);
	save_volume4D(data, string(outputFile));

	// Baseline Subtraction
	for (int t = 0; t < data.tsize(); t++)
	{
		if (interval.IsActivationBlockStart(t + 1))
		{
			vol = data[t];
			withBaseline = 1;
		}

		// actual subtraction
		if (withBaseline) data[t] = data[t] - vol;
		else data[t] = 0;
	}

	// gaussian smoothing
	float xdim = data.xdim();
	float ydim = data.ydim();
	float zdim = data.zdim();

	volume<float> kernel(box_kernel(3, 3, 3));
	kernel = gaussian_kernel3D(FWHM, xdim, ydim, zdim);
	generic_convolve(data, kernel, true, true);
	sprintf(outputFile, "%s\\%s.nii", dirSuj, dirRun);
	save_volume4D(data, string(outputFile));
}

void processDirectory(char *dirStudy, char *subject, char *feedback, float FWHM)
{
	volume4D<float> data;
	int nvols = 296;

	volume<float> vol;
	char fname[1024];
	char dirRun[1024];
	char dirSuj[1024];

	sprintf(dirSuj, "%s\\MNI\\%s", dirStudy, subject);
	for (int r = 1; r < 5; r++)
	{
		char regressFile[1024], outName[1024], designFile[1024];

		sprintf(dirRun, "RUN%.2dNFB0%c", r, feedback[r]);

		sprintf(outName, "%s\\%s.nii", dirSuj, dirRun);
		remove(outName);
		sprintf(designFile, "%s\\DESENHO_NFB0%c.txt", dirStudy, feedback[r]);
		sprintf(regressFile, "%s\\GLM\\%s_MP_regress.mat", dirSuj, dirRun);

		DesignObject interval;
		interval.readDesignFile(designFile);
		strcpy(interval.baselineCondition, "NEUTRO");

		if (1)
		{
			for (int t = 1; t <= nvols; t++)
			{
				sprintf(fname, "%s\\%s\\DRIN-%.5d.nii", dirSuj, dirRun, t);
				read_volume(vol, string(fname));
				if (t == 1)
				{
					data.reinitialize(vol.xsize(), vol.ysize(), vol.zsize(), nvols);
					data.copyproperties(vol);
				}
				data[t - 1] = vol;
			}
			processData(data, FWHM, regressFile, interval, dirSuj, dirRun);
			runMelodic(dirSuj, dirRun);
		}
		else read_volume4D(data, string(outName));
		//svmTrain(data, interval, dirSuj, dirRun);
	}
}

int maskSize(NEWIMAGE::volume<float> &mask)
{
	int i = 0;
	for (int z = 0; z < mask.zsize(); z++)
		for (int y = 0; y < mask.ysize(); y++)
			for (int x = 0; x < mask.xsize(); x++)
				if (mask.value(x, y, z) > 0)
				{
					i++;
				}
	return i;
}

int betCheck(volume<float> &volume)
{
	int nvox = volume.nvoxels();
	int nvoxnz = maskSize(volume);
	if ((float)(nvoxnz / nvox) > 0.9) return 0;
	else return 1;

}

void clean_up_edge(volume<float> &vol, volume<float> &mask, char *outputImg)
{
	//# does some despiking filtering to clean up the edge of the fieldmap
	//# args are : <fmap> <mask> <tmpnam>
	//outfile = $1
	//maskim = $2
	//tmpnm = $3

	mkdir("workdir");
	char CmdLn[1024];
	// fugue --loadfmap=${outfile} --savefmap=${tmpnm}_tmp_fmapfilt --mask=${maskim} --despike --despikethreshold=2.1
	save_volume(vol, "workdir/tmp_fmap");
	save_volume(vol, string(outputImg));
	sprintf(CmdLn, "fugue --loadfmap=workdir/tmp_fmap --savefmap=workdir/tmp_fmapfilt --mask=workdir/mask --despike --despikethreshold=2.1");
	fugue(CmdLn);

	// fslmaths ${maskim} -kernel 2D -ero ${tmpnm}_tmp_eromask
	sprintf(CmdLn, "fslmaths workdir/mask -kernel 2D -ero workdir/tmp_eromask");
	fslmaths(CmdLn);

	// fslmaths ${maskim} -sub ${tmpnm}_tmp_eromask -thr 0.5 -bin ${tmpnm}_tmp_edgemask
	sprintf(CmdLn, "fslmaths workdir/mask -sub workdir/tmp_eromask -thr 0.5 -bin workdir/tmp_edgemask");
	fslmaths(CmdLn);

	// fslmaths ${tmpnm}_tmp_fmapfilt -mas ${tmpnm}_tmp_edgemask ${tmpnm}_tmp_fmapfiltedge
	sprintf(CmdLn, "fslmaths workdir/tmp_fmapfilt -mas workdir/tmp_edgemask workdir/tmp_fmapfiltedge");
	fslmaths(CmdLn);

	// fslmaths ${outfile} -mas ${tmpnm}_tmp_eromask -add ${tmpnm}_tmp_fmapfiltedge ${outfile}
	sprintf(CmdLn, "fslmaths %s -mas workdir/tmp_eromask -add workdir/tmp_fmapfiltedge %s", outputImg, outputImg);
	fslmaths(CmdLn);
}


void demean_image(volume<float> &vol, volume<float> &mask) 
{
	volume<float> tempMask;
	tempMask = mask;
	tempMask.binarise(0, tempMask.max() + 1, exclusive);

	// Masking volume
	vol *= tempMask;

	// Getting the mean value ( 50 percentile )
	vol -= vol.percentile(0.5, vol);
}

void fsl_prepare_fieldmap(char *phaseImage, char *magImage, char *outputImg, float deltaTE)
{
	char CmdLn[1024];
	if ((deltaTE < 0.1) || (deltaTE > 10.0)) return;
	volume<float> phaseVol, magVol;

	read_volume(phaseVol, string(phaseImage));
	read_volume(magVol, string(magImage));

	// dimensions check
	if ((phaseVol.xdim() != magVol.xdim()) || (phaseVol.ydim() != magVol.ydim()) || (phaseVol.zdim() != magVol.zdim())
		|| (phaseVol.xsize() != magVol.xsize()) || (phaseVol.ysize() != magVol.ysize()) || (phaseVol.zsize() != magVol.zsize())) return;

	// range check
	float minIntensid = phaseVol.min();
	float maxIntensid = phaseVol.max();
	float range = maxIntensid - minIntensid;
	float nrange = range / 4096.0;
	if (nrange < 2.1)
	{
		if (nrange > 1.9)
		{
			phaseVol /= 2;
		}
	}
	// bet check
	if (!betCheck(magVol)) return;
	mkdir("workdir");


	// Make brain mask
	volume<float> mask;
	mask = magVol;
	mask.threshold(0.00000001);
	mask.binarise(0, mask.max() + 1, exclusive);

	// Convert phasemap to radians
	// fslmaths ${newphaseroot} -div 2048 -sub 1 -mul 3.14159 -mas ${maskim} ${tmpnm}_tmp_ph_radians -odt float
	save_volume(phaseVol, "workdir/radians");
	save_volume(magVol, "workdir/mag");
	save_volume(mask, "workdir/mask");
	sprintf(CmdLn, "fslmaths workdir/radians -div 2048 -sub 1 -mul 3.14159 -mas workdir/mask workdir/ph_radians -odt float");
	fslmaths(CmdLn);

	//# Unwrap phasemap
	// prelude - p ${ tmpnm }_tmp_ph_radians - a ${ absroot } -m ${ maskim } -o ${ tmpnm }_tmp_ph_radians_unwrapped - v
	sprintf(CmdLn, "prelude -p workdir/ph_radians -a workdir/mag -m workdir/mask -o workdir/unwrapped -v");
	prelude(CmdLn);

	// Convert to rads / sec(dTE is echo time difference)
	// asym = `echo $dTE / 1000 | bc - l`
	// fslmaths ${tmpnm}_tmp_ph_radians_unwrapped -div $asym ${tmpnm}_tmp_ph_rps -odt float
	sprintf(CmdLn, "fslmaths workdir/unwrapped -div %f workdir/tmp_ph_rps -odt float", (deltaTE / 1000));
	fslmaths(CmdLn);

	// Call FUGUE to extrapolate from mask(fill holes, etc)
	// fugue --loadfmap=${ tmpnm }_tmp_ph_rps --mask=${maskim} --savefmap=$outfile
	sprintf(CmdLn, "fugue --loadfmap=workdir/tmp_ph_rps --mask=workdir/mask --savefmap=%s", outputImg);
	fugue(CmdLn);

	//# Demean to avoid gross shifting
	// demean_image ${outfile} ${maskim} ${tmpnm}
	read_volume(phaseVol, string(outputImg));
	demean_image(phaseVol, mask);

	//# Clean up edge voxels
	// clean_up_edge ${outfile} ${maskim} ${tmpnm}
	clean_up_edge(phaseVol, mask, outputImg);
}

void getFieldMapInfo(char *fileName, float &dwellTime, int &direction, float &TE)
{
	FSLIO *OP = FslOpen(fileName, "r");
	if (OP != NULL)
	{
		char description[100];
		strcpy(description, OP->niftiptr->descrip);

		strings tags, auxtags;
		splitLine(string(description), tags, ';');
		if (tags.size() > 0)
		{
			splitLine(tags[0], auxtags, '=');
			TE = stod(auxtags[1]);
		}
		if (tags.size() > 2)
		{
			splitLine(tags[2], auxtags, '=');
			direction = stoi(auxtags[1]);
		}
		if (tags.size() > 3)
		{
			splitLine(tags[3], auxtags, '=');
			dwellTime = stod(auxtags[1]) / 1000;
		}
		FslClose(OP);
		FslFree(OP);
	}
}

void brainThreshold(char *functional, char *output, float brainThres)
{
	float p2, p98;
	float thres;
	char cmd[500];
	char quote[10];
	volume4D<float> Vol;

	sprintf(quote, "\"");

	// 2th percentil and 98th percentil calculation
	read_volume4D(Vol, string(functional));
	p2 = Vol.percentile(0.02);
	p98 = Vol.percentile(0.98);

	// threshold calculation. Minimum value equals a percentage from the difference between the 98th and 2th percentile. Note the brain threshold parameter is in the range [1..100]
	thres = p2 + (brainThres / 100) * (p98 - p2);

	// applying the threshold and get the minimum value per voxel, in case of a more than one volume in file.
	sprintf(cmd, "fslmaths %s%s%s -thr %f -Tmin -bin %s%s%s -odt char", quote, functional, quote, thres, quote, output, quote);
	fslmaths(cmd);
}

void prepareFieldMap(char *reference, char *phaseImage, char *magImage, char *fieldMap, float deltaTE)
{
	// bet functional
	char CmdLn[1024];
	mkdir("workdir");
	
	sprintf(CmdLn, "bet %s workdir/%s -f 0.3", reference, "reference_bet");
	bet(CmdLn);
	sprintf(CmdLn, "bet workdir/%s workdir/%s -f 0.3", "reference_bet", "reference_bet");
	bet(CmdLn);

	// bet magImage
	sprintf(CmdLn, "bet %s workdir/%s -f 0.3", magImage, "magimage_bet");
	bet(CmdLn);
	sprintf(CmdLn, "bet workdir/%s workdir/%s -f 0.3", "magimage_bet", "magimage_bet");
	bet(CmdLn);

	
	sprintf(CmdLn, "flirt -ref workdir/%s -in workdir/%s -dof 7 -out workdir/%s", "magimage_bet", "reference_bet", "magimage_mask");
	flirt(CmdLn);

	char magnitudeMask[1024], magnitudeThresholded[1024], magnitudeBet[1024], magnitudeMaskFinal[1024];

	sprintf(magnitudeMask, "workdir/%s", "magimage_mask");
	sprintf(magnitudeThresholded, "workdir/%s", "magimage_thres");
	sprintf(magnitudeBet, "workdir/%s", "magimage_bet");
	sprintf(magnitudeMaskFinal, "workdir/%s", "magimage_mask_final");

	brainThreshold(magnitudeMask, magnitudeThresholded, 15);

	sprintf(CmdLn, "fslmaths %s -mas %s %s", magnitudeBet, magnitudeThresholded, magnitudeMaskFinal);
	fslmaths(CmdLn);
	fsl_prepare_fieldmap(phaseImage, magnitudeMaskFinal, fieldMap, deltaTE);

	//removeDirectory("workdir");
}

void applyFieldMap(char *fileName, char *fieldMap, char *output, float dwellTime, int direction)
{
	char CmdLn[2048], unwarpdir[10];

	if (direction == 1) sprintf(unwarpdir, "y-");
	else sprintf(unwarpdir, "y");
    sprintf(CmdLn, "fugue -i %s --dwell=%f --loadfmap=%s -u %s --unwarpdir=%s", fileName, dwellTime, fieldMap, output, unwarpdir);
	fugue(CmdLn);
}

void applyFieldMap(char *fileName, char *fieldMap, char *output)
{
	float dwellTime;
	int direction;
	float TE;
	getFieldMapInfo(fileName, dwellTime, direction, TE);
	applyFieldMap(fileName, fieldMap, output, dwellTime, direction);
}
