#include "fMRIProcessing.h"

int createSmallSliceReport(char *fileName, char *rootDir = NULL)
{
	char cmd[2040], dirPai[2048], tempFile[2048];

	if (!rootDir) extractFilePath(fileName, dirPai);
	else strcpy(dirPai, rootDir);

	extractFileName(fileName, tempFile);
	string inputName = rootDir + string("\\") + tempFile;
	make_basename(inputName);

	if (chdir(dirPai) == 0)
	{
		sprintf(cmd, "slicer %s -s 3 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png", fileName);
		slicer(cmd);

		sprintf(cmd, "pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png %s.png", inputName.c_str());
		pngappend(cmd);
		// apagar todos os pngs
	}
	return 0;
}

int createSliceReport(char *file1, char *file2, char *pngFile, char *rootDir = NULL)
{
	char cmd[2040], dirPai[2048], tempFile[2048];

	if (!rootDir) extractFilePath(file1, dirPai);
	else strcpy(dirPai, rootDir);

	extractFileName(file1, tempFile);
	string inputName = rootDir + string("\\") + tempFile;
	make_basename(inputName);

	if (chdir(dirPai) == 0)
	{
		sprintf(cmd, "slicer %s %s -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png", file1, file2);
		slicer(cmd);

		sprintf(cmd, "pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png %s_1.png", inputName.c_str());
		pngappend(cmd);

		sprintf(cmd, "slicer %s %s -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png", file2, file1);
		slicer(cmd);

		sprintf(cmd, "pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png %s_2.png", inputName.c_str());
		pngappend(cmd);

		sprintf(cmd, "pngappend %s_1.png - %s_2.png %s.png", inputName.c_str(), inputName.c_str(), inputName.c_str());
		pngappend(cmd);

		sprintf(pngFile, "%s.png", inputName.c_str());
	}

	// apagar todos os pngs
	return 0;
}

int applyTransformation(char *volume, char *output, char *reference, char *matrixFile, int nn)
{
	char cmd[2040];
	char interpMode[20];

	if (nn) strcpy(interpMode, "nearestneighbour");
	else strcpy(interpMode, "trilinear");

	//	sprintf(cmd, "flirt -in %s -ref %s -out %s -applyxfm -paddingsize 10 -init %s -interp %s -verbose 5", volume, reference, output, matrixFile, interpMode);
	sprintf(cmd, "flirt -in %s -ref %s -out %s -applyxfm -paddingsize 10 -init %s -interp %s", volume, reference, output, matrixFile, interpMode);
	flirt(cmd);

	return 0;
}

int mniRegistration(char *volume, char *highres, char *reference, char *pngFile, char *matrixFile)
{
	char cmd[2040], dirPai[2048], tempFile[2048], file1[2048], file2[2048];
	char dof[] = "12";
	char interpMode[] = "trilinear";
	float angle = 90;

	extractFilePath(volume, dirPai);

	extractFileName(volume, file1);
	extractFileName(reference, file2);

	string file1S = file1;
	make_basename(file1S);

	string file2S = file2;
	make_basename(file2S);

	string outputName = file1S + "2" + file2S;

	char regDir[BUFF_SIZE], fileName[BUFF_SIZE];

	sprintf(regDir, "%s\\reg", dirPai);
	mkdir(regDir);

	outputName = string(regDir) + "\\" + outputName;

	sprintf(cmd, "flirt -in %s -ref %s -out %s -omat %s.mat -cost corratio -dof %s -searchrx -%f %f -searchry -%f %f -searchrz -%f %f -interp %s", volume, reference, outputName.c_str(), outputName.c_str(), dof, angle, angle, angle, angle, angle, angle, interpMode);
	flirt(cmd);

	sprintf(matrixFile, "%s.mat", outputName.c_str());
	createSmallSliceReport(reference, regDir);
	createSliceReport((char *)outputName.c_str(), reference, pngFile, regDir);

	return 0;
}
