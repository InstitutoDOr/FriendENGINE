#include "fMRIProcessing.h"
#include "intervals.h"
#include "vectorUtils.h"
#include "fMRIUtils.h"

int regressMovementsParams4D(char *intervalFile, char *condBASAL, char *iname, char *oname, int runSize, volume4D<float> &data, int includeMotionParameters = 1)
{
	char auxString[2048], glmDir[2048];
	string confoundFile, fsfFile;

	confoundFile = iname;
	make_basename(confoundFile);
	confoundFile = confoundFile + ".par";

	stringstream CmdLn;
	DesignObject interval;
	strcpy(interval.baselineCondition, condBASAL);
	interval.readDesignFile(intervalFile);

	extractFilePath(iname, auxString);
	sprintf(glmDir, "%s\\glm", auxString);
	mkdir(glmDir);

	fsfFile = glmDir;
	fsfFile = fsfFile + "\\design";

	// generating the design matrix file
	if (interval.conditionNames.size() > 0)
	{
		fprintf(stderr, "Generating conditions files.\n");
		interval.generateConditionsBoxCar(glmDir);
		fprintf(stderr, "Generating FSF file.\n");

		strcpy(interval.glmDir, glmDir);
		interval.generateFSFFile((char *)(fsfFile + ".fsf").c_str(), runSize, includeMotionParameters, 0);

		fprintf(stderr, "Running feat_model.\n");
		CmdLn.str("");

		CmdLn << "feat_model " << fsfFile;
		if (includeMotionParameters) CmdLn << " " << confoundFile;
		chdir(glmDir);

		//fprintf(stderr, "%s\n", (char *)CmdLn.str().c_str());
		feat_model((char *)CmdLn.str().c_str());
	}
	regressOutDesign((char *)(fsfFile + ".mat").c_str(), (char *)(fsfFile + "_regressed.txt").c_str(), interval.conditionNames.size() - 1);
	regressOutData(data, (char *)(fsfFile + "_regressed.txt").c_str());
	return 1;
}

int motionCorrection4D(char *volume4D, char *reference, char *regressorFile, char *mcsuffix, char *filesuffix)
{
	char cmd[2040], dirPai[2048], tempFile[2048];
	extractFilePath(volume4D, dirPai);
	string inputfname = volume4D;
	make_basename(inputfname);

	inputfname = inputfname + "_mcf";
	if (filesuffix)
		inputfname = inputfname + "_" + string(filesuffix);

	char mcDir[BUFF_SIZE], fileName[BUFF_SIZE];

	extractFileName((char *)inputfname.c_str(), fileName);

	if (mcsuffix) sprintf(mcDir, "%s\\mc_%s", dirPai, mcsuffix);
	else sprintf(mcDir, "%s\\mc", dirPai);

	sprintf(tempFile, "%s\\%s", mcDir, fileName);

	// criar pasta mc
	mkdir(mcDir);

	fprintf(stderr, "Motion Correction\n");

	sprintf(cmd, "mcflirt -in %s -out %s -mats -plots -reffile %s -rmsrel -rmsabs -spline_final", volume4D, tempFile, reference);
	mcflirt(cmd);
	sprintf(regressorFile, "%s.par", tempFile);

	// criar trans.png
	if (filesuffix) sprintf(cmd, "fsl_tsplot -i %s.par -t \"MCFLIRT estimated translations (mm)\" -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o %s\\trans_%s.png", tempFile, mcDir, filesuffix);
	else sprintf(cmd, "fsl_tsplot -i %s.par -t \"MCFLIRT estimated translations (mm)\" -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o %s\\trans.png", tempFile, mcDir);
	fprintf(stderr, "%s\n", cmd);
	fsl_tsplot(cmd);

	// criar rot.png
	if (filesuffix) sprintf(cmd, "fsl_tsplot -i %s.par -t \"MCFLIRT estimated rotations (radians)\" -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o %s\\rot_%s.png", tempFile, mcDir, filesuffix);
	else sprintf(cmd, "fsl_tsplot -i %s.par -t \"MCFLIRT estimated rotations (radians)\" -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o %s\\rot.png", tempFile, mcDir);
	fprintf(stderr, "%s\n", cmd);
	fsl_tsplot(cmd);

	// disp.png rms
	if (filesuffix) sprintf(cmd, "fsl_tsplot -i %s_abs.rms,%s_rel.rms -t \"MCFLIRT estimated mean displacement (mm)\" -u 1 --start=0 --finish=0 -w 640 -h 144 -a absolute,relative -o %s\\disp_%s.png", tempFile, tempFile, mcDir, filesuffix);
	else sprintf(cmd, "fsl_tsplot -i %s_abs.rms,%s_rel.rms -t \"MCFLIRT estimated mean displacement (mm)\" -u 1 --start=0 --finish=0 -w 640 -h 144 -a absolute,relative -o %s\\disp.png", tempFile, tempFile, mcDir);
	fprintf(stderr, "%s\n", cmd);
	fsl_tsplot(cmd);

	sprintf(cmd, "fslmaths %s %s", tempFile, inputfname.c_str());
	fslmaths(cmd);
	fileExists(tempFile);
	remove(tempFile);
	return 0;
}

void create4Dvolume(char *templateFile, char *output)
{
	int numVolumes = 1;
	int saida = 0, lido = 0;
	char fileName[2048];

	while (!saida)
	{
		sprintf(fileName, "%s%d.nii", templateFile, numVolumes);
		if (fileExists(fileName)) numVolumes++;
		else saida = 1;
	}

	volume4D<float> thr4D;
	for (int t = 1; t < numVolumes; t++)
	{
		sprintf(fileName, "%s%d.nii", templateFile, t);
		volume<float> vol;
		if (fileExists(fileName))
		{
			read_volume(vol, string(fileName));
			if (t == 1)
			{
				thr4D.copyproperties(vol);
				thr4D.reinitialize(vol.xsize(), vol.ysize(), vol.zsize(), numVolumes - 1);
				lido = 1;
			}
			thr4D[t - 1] = vol;
		}
	}
	if (lido)
		save_volume4D(thr4D, string(output));
}

////////////////////////////////// 4D Volume processing /////////////////////////////////////////////////////////////////////////////
int processSubject(char *projectdir, char *subject, char *standard, vector<string>&runs, int window)
{
	int debug = 0;
	char subjectdir[2048], mapsdir[2048], qcdir[2048];
	sprintf(subjectdir, "%s\\%s", projectdir, subject);
	mkdir(subjectdir);

	sprintf(mapsdir, "%s\\Maps", projectdir);
	mkdir(mapsdir);

	sprintf(qcdir, "%s\\QC", projectdir);
	mkdir(qcdir);

	if (runs.size())
	{
		char toFile[2048], fromFile[2048], matrixFile[2048];

		// criar a referencia
		char reference[2048];
		sprintf(reference, "%s\\reference", subjectdir);
		extractVolume((char *)runs[0].c_str(), reference, 0, getDataType((char *)runs[0].c_str()));

		// referencia mni
		mniRegistration(reference, "", standard, fromFile, matrixFile);

		sprintf(toFile, "%s\\%s_registration.png", qcdir, subject);
		copyFile(fromFile, toFile);
		for (int r = 1; r <= runs.size(); r++)
		{
			char runDir[2048], runFile[2048], cmd[2048];
			char actualFile[2048];

			if (!debug)
			{
				char regressorFile[2048];
				// criar um diretorio para cada run
				sprintf(runDir, "%s\\RUN%.2d", subjectdir, r);
				mkdir(runDir);

				sprintf(runFile, "%s\\RUN%.2d", runDir, r);

				// copiando o run
				sprintf(cmd, "fslmaths %s %s", runs[r - 1].c_str(), runFile);
				fslmaths(cmd);

				// fazer o motion correction
				motionCorrection4D(runFile, reference, regressorFile);

				volume4D<float> data;
				sprintf(actualFile, "%s_mcf", runFile);
				read_volume4D(data, string(actualFile));

				// regressar os parametros de movimento
				regressOutData(data, regressorFile);

				// Sliding window
				if (window > 0)
				{
					// Sliding Window
					fprintf(stderr, "Sliding Window\n");
					estimateActivation4D(data, window);

					sprintf(actualFile, "%s_mcfs", runFile);
					save_volume4D(data, actualFile);
				}
				runMelodic4D(actualFile);
			}
			if (debug)
				sprintf(actualFile, "E:\\projetos\\PRJ1016_PRIMING\\SUBJ045\\RUN%.2d\\RUN%.2d_mcfs", r, r);

			sprintf(fromFile, "%s.ica\\filtered_func_data.nii.ica\\melodic_IC.nii.gz", actualFile);
			sprintf(toFile, "%s\\%s_RUN%.2d_IC_maps.nii.gz", mapsdir, subject, r);
			copyFile(fromFile, toFile);

			sprintf(fromFile, "%s.ica\\filtered_func_data.nii.ica\\melodic_mix", actualFile);
			sprintf(toFile, "%s\\%s_RUN%.2d_IC_timecourses.txt", mapsdir, subject, r);
			copyFile(fromFile, toFile);

			sprintf(fromFile, "%s.ica\\filtered_func_data.nii.ica\\stats\\thresh_zstat", actualFile);
			sprintf(toFile, "%s\\%s_RUN%.2d_IC_thr_maps.nii.gz", mapsdir, subject, r);
			create4Dvolume(fromFile, toFile);
			sprintf(fromFile, "%s", toFile);
			sprintf(toFile, "%s\\%s_RUN%.2d_IC_thr_maps_MNI.nii.gz", mapsdir, subject, r);
			applyTransformation(fromFile, toFile, standard, matrixFile);
		}
	}
	return 0;
}

int processVolumeRun4D(char *fileName, char *intervalFile, char* condBASAL, int window = 0)
{
	char cmd[2040], dirPai[2048];
	char tempString[2048];

	int runSize;
	short dataType;

	volumeProperties(fileName, runSize, dataType);
	// MOTION CORRECTION
	if (1)
	{
		fprintf(stderr, "Motion Correction\n");
		extractFilePath(fileName, dirPai);
		sprintf(tempString, "%s\\example_func", dirPai);
		int middleIdx = 0;
		extractVolume(fileName, tempString, middleIdx, dataType);
		sprintf(cmd, "mcflirt -in %s -mats -plots -reffile %s -rmsrel -rmsabs -spline_final", fileName, tempString);
		mcflirt(cmd);
	}
	string inputfname = fileName;
	string outputfname;

	make_basename(inputfname);
	outputfname = inputfname + "_mcf";

	volume4D<float> data;

	read_volume4D(data, outputfname);

	if (fileExists(intervalFile))
	{
		// Regress Motion Params
		fprintf(stderr, "Regress Out Movement Params\n");
		regressMovementsParams4D(intervalFile, condBASAL, (char *)outputfname.c_str(), (char *)(inputfname + "_mcfr").c_str(), runSize, data);
		save_volume4D(data, inputfname + "_mcfr");
	}

	if (window > 0)
	{
		// Sliding Window
		fprintf(stderr, "Sliding Window\n");
		estimateActivation4D(data, window);
	}
	save_volume4D(data, inputfname + "_mcfrs");
	runMelodic4D((char *)(inputfname + "_mcfrs").c_str());
	return 0;
}

void detrendSG(char *filename, char *saida, int w, int degree)
{
	volume4D<float> volume;
	read_volume4D(volume, string(filename));

	float_vect unfiltered, filtered;
	unfiltered.resize(volume.tsize());

	for (int z = 0; z < volume.zsize(); z++)
	{
		for (int y = 0; y < volume.ysize(); y++)
		{
			for (int x = 0; x < volume.xsize(); x++)
			{
				for (int t = 0; t < volume.tsize(); t++) unfiltered[t] = volume.value(x, y, z, t);
				filtered = detrendVector(unfiltered, w, degree);
				for (int t = 0; t < volume.tsize(); t++) volume.value(x, y, z, t) = filtered[t];
			}
		}
	}
	save_volume4D(volume, saida);
}
