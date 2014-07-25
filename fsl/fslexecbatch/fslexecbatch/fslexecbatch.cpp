// fslexecbatch.cpp : Defines the entry point for the console application.
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	char cmd[500];
	char dir[500];
	char dirresults[500];
	int offset;
	double threshold;
	char *sufixos[] = {"gito\\fig\\", "gito\\imag\\", "ricardo\\fig\\", "ricardo\\imag\\", "JEFFERSON\\fig\\", "JEFFERSON\\imag\\"};

	for (int i=3;i<6;i++)
	{
	   strcpy(dirresults, "C:\\Instituto\\x64\\Release\\resultados\\");
	   strcat(dirresults, sufixos[i]);

	   strcpy(dir, "C:\\projetos\\");
	   strcat(dir, sufixos[i]);
	   for (threshold=2;threshold<7;threshold=threshold+0.5)
	   {
		  sprintf(cmd, "c:\\Instituto\\x64\\Release\\fslexec unemascaras %sglmsaidat.nii %sfeatures.nii %.2f", dir, dir, threshold);
		  system(cmd);
		  for (offset = 2;offset < 6; offset++)
		  {
			 sprintf(cmd, "c:\\Instituto\\x64\\Release\\fslexec svm %sbet_gauss_box.nii %sdesenho.txt %sfeatures.nii %d %scmatrix_THRESHOLD%.2f_OFFSET%d.txt", dir, dir, dir, offset, dirresults, threshold, offset);
			 system(cmd);
		  }
	   }
	}
	return 0;
}

