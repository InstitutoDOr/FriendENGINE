#include "broccoliInterface.h"

using namespace NEWIMAGE;

int sameBBox(char *v1, char *v2)
{
	volume<float>T1v, MNIv;

	read_volume(T1v, string(v1));
	read_volume(MNIv, string(v2));

	int x1 = T1v.maxx() - T1v.minx() + 1;
	int y1 = T1v.maxy() - T1v.miny() + 1;
	int z1 = T1v.maxz() - T1v.minz() + 1;

	int x2 = MNIv.maxx() - MNIv.minx() + 1;
	int y2 = MNIv.maxy() - MNIv.miny() + 1;
	int z2 = MNIv.maxz() - MNIv.minz() + 1;

	int newx = max(x1, x2);
	int newy = max(y1, y2);
	int newz = max(z1, z2);

	char n1[500], n2[500];

	volume<float>T1o(newx, newy, newz);
	pad(T1v, T1o);
	changeFileExt(v1, "_bb.nii", n1);
	save_volume(T1o, string(n1));

	volume<float>MNIo(newx, newy, newz);
	pad(MNIv, MNIo);
	changeFileExt(v2, "_bb.nii", n2);
	save_volume(MNIo, string(n2));

	return 0;
}

void xCut(char *inputFilename, char *outputFilename, int axis, int ini, int end)
{
	volume<float>inputVol;

	read_volume(inputVol, string(inputFilename));

	volume<float>outputVol(inputVol.xsize(), inputVol.ysize(), inputVol.zsize());
	outputVol.copyproperties(inputVol);

	int inix = inputVol.minx();
	int iniy = inputVol.miny();
	int iniz = inputVol.minz();

	int endx = inputVol.maxx();
	int endy = inputVol.maxy();
	int endz = inputVol.maxz();

	if (axis == 1)
	{
		if (ini > inix) inix = ini;
		if ((end) && (end < endx)) endx = end;
	}
	else if (axis == 2)
	{
		if (ini > iniy) iniy = ini;
		if ((end) && (end < endy)) endy = end;
	}
	else if (axis == 3)
	{
		if (ini > iniz) iniz = ini;
		if ((end) && (end < endz)) endz = end;
	}

	outputVol = 0;
	for (int z = iniz; z <= endz; z++)
		for (int y = iniy; y <= endy; y++)
			for (int x = inix; x <= endx; x++)	outputVol(x, y, z) = inputVol(x, y, z);

	save_volume(outputVol, string(outputFilename));
}

int main(int argc, char **argv)
{
	char arqAxial[500], arqAxialBB[500], arqAxialBBBet[500], T1[500], Template[500], MNIn[500];
	strcpy(MNIn, "E:\\FriendEngine\\BROCCOLIENGINE\\Engine\\Release\\MNI152_T1_1mm_brain.nii");
	strcpy(T1, "E:\\FriendEngine\\BROCCOLIENGINE\\Engine\\Release\\RAI.nii");
	strcpy(arqAxial, "E:\\FriendEngine\\BROCCOLIENGINE\\Engine\\Release\\RAI_ax.nii");
	strcpy(arqAxialBB, "E:\\FriendEngine\\BROCCOLIENGINE\\Engine\\Release\\RAI_ax_bb.nii");
	strcpy(arqAxialBBBet, "E:\\FriendEngine\\BROCCOLIENGINE\\Engine\\Release\\RAI_ax_bb_bet.nii");

	if (0)
	{
		resampleVolume(T1, arqAxial, 1, 1, 1, 2, 0);
		axial(arqAxial, arqAxial);
		xCut(arqAxial, arqAxialBB, 3, 70, 0);

		char cmd[1024];
		sprintf(cmd, "bet %s %s -f 0.3", arqAxialBB, arqAxialBBBet);
		bet(cmd);

//		sprintf(cmd, "bet %s %s -f 0.3", arqAxialBBBet, arqAxialBBBet);
//		bet(cmd);
		BROCCOLILinearRegistration(arqAxialBBBet, "E:\\FriendEngine\\BROCCOLIENGINE\\Engine\\Release\\MNI152_T1_1mm_brain.nii", "E:\\FriendEngine\\BROCCOLIENGINE\\Engine\\Release\\saida.nii");
	}

	if (0)
	{
		double startTime = GetWallTime();
		BROCCOLIMotionCorrection("E:\\BTeste\\RFI.nii", "E:\\BTeste\\RUN01\\FUNC4D_", "E:\\BTeste\\confounds.txt", 304);
		double endTime = GetWallTime();

		printf("\nIt took %f seconds to run the motion correction\n", (float)(endTime - startTime));
	}
	else
	{
		char outputFile[500], prefix[500] = "_RFI2", prefix2[500], prefix3[500];
		char cmd[1024];

		if (0)
		{
			sprintf(outputFile, "E:\\BTeste\\saida_RFI2.nii");
			MniToSubject("E:\\BTeste\\RFI.nii", "E:\\projetos\\NFB_VR\\Mask_VR_ROI.nii", "E:\\BTeste\\MNI152_T1_1mm_brain.nii", outputFile, prefix);

			for (int t = 1; t <= 304; t++)
			{
				sprintf(cmd, "flirt -ref E:\\BTeste\\MNI152_T1_1mm_brain.nii -in E:\\BTeste\\RUN01\\FUNC4D_%.5d.nii -o E:\\BTeste\\RUN01MNI\\FUNC4D_%.5d.nii -applyxfm -init E:\\BTeste\\RFI2MNI_RFI2.mat -interp nearestneighbour", t, t);
				flirt(cmd);
			}
			return 0;
		}

		BROCCOLIEngineInterface interface;
		interface.createBROCCOLIObject(0, 0);

		if (0)
		{
			strcpy(prefix, "E:\\BTeste\\RUN01MNI\\FUNC4D_");
			strcpy(prefix2, "E:\\BTeste\\RUN01.1\\FUNC4D_");
			strcpy(prefix3, "E:\\BTeste\\RUN01.2\\FUNC4D_");
			interface.setRFI("E:\\BTeste\\RFI2MNI_RFI2.nii");
		}
		else
		{
			strcpy(prefix, "E:\\BTeste\\RUN01\\FUNC4D_");
			strcpy(prefix2, "E:\\BTeste\\RUN01.1\\FUNC4D_");
			strcpy(prefix3, "E:\\BTeste\\RUN01.2\\FUNC4D_");
			interface.setRFI("E:\\BTeste\\RFI.nii");
		}
		interface.preparePipeline(5);

		//// BROCCOLI PIPELINE //////////////////////////////////
		double startTime = GetWallTime();
		for (int t = 1; t <= 304; t++)
		{
			char filename[500], outputName[500];

			sprintf(filename, "%s%.5d.nii", prefix, t);
			sprintf(outputName, "%s%.5d_mc_g_1.nii", prefix2, t);
			interface.pipelineEngine(filename, outputName);
		}
		double endTime = GetWallTime();
		interface.deallocateBROCCOLIObject();

		printf("\nIt took %f seconds to run the engine pipeline.\n", (float)(endTime - startTime));

		if (1)
		{
			//// NORMAL PIPELINE //////////////////////////////////
			startTime = GetWallTime();
			for (int t = 1; t <= 304; t++)
			{
				char filename[500], outputName[500];

				sprintf(filename, "%s%.5d.nii", prefix, t);
				sprintf(outputName, "%s%.5d_mc_g_2.nii", prefix3, t);
				interface.pipelineNormalEngine(filename, outputName);
			}
			endTime = GetWallTime();
			printf("\nIt took %f seconds to run the normal engine pipeline.\n", (float)(endTime - startTime));
		}
	}
}