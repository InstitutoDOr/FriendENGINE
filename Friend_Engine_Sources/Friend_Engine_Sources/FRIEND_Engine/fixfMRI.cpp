#include "fMRIProcessing.h"

void createEdges(char *melodicDir)
{
	char fixDir[1024];

	sprintf(fixDir, "%s\\fix", melodicDir);
	mkdir(fixDir);

	// create edge masks
	char cmd[2048];

	sprintf(cmd, "fslmaths %s\\mask -eroF %s\\maske1", melodicDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\maske1 -eroF %s\\maske2", fixDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\maske2 -eroF %s\\maske3", fixDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\maske3 -eroF %s\\maske4", fixDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\maske4 -eroF %s\\maske5", fixDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\mask -sub %s\\maske1 %s\\edge1", melodicDir, fixDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\mask -sub %s\\maske2 %s\\edge2", melodicDir, fixDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\mask -sub %s\\maske3 %s\\edge3", melodicDir, fixDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\mask -sub %s\\maske4 %s\\edge4", melodicDir, fixDir, fixDir);
	fslmaths(cmd);

	sprintf(cmd, "fslmaths %s\\mask -sub %s\\maske5 %s\\edge5", melodicDir, fixDir, fixDir);
	fslmaths(cmd);
}

void applyFast(char *melodicDir)
{
	char fixDir[1024], cmd[2048];
	char priorTemplate[] = "H:\\FriendEngine\\SIC\\x64\\Release\\tissuepriors\\avg152T1";

	sprintf(fixDir, "%s\\fix", melodicDir);
	mkdir(fixDir);
	sprintf(cmd, "fast -t 1 -o %s\\fastsg -A %s_csf %s_gray %s_white %s\\reg\\highres ", fixDir, priorTemplate, priorTemplate, priorTemplate, melodicDir);
	fast(cmd);

	//sprintf(cmd, "flirt -in %s\\fastsg_pveseg -ref %s\\reg\\reference -applyxfm -init %s\\reg\\highres2example_func.mat -interp nearestneighbour -out %s\\hr2exf", fixDir, melodicDir, melodicDir, fixDir);
	//flirt(cmd);
}


void fix(char *melodicDir)
{
	char fixDir[1024];

	sprintf(fixDir, "%s\\fix", melodicDir);
	mkdir(fixDir);

	// create edge masks
	createEdges(melodicDir);

	// apply Fast
	applyFast(melodicDir);
}
