/*
#include <stdio.h>
#include "engine.h"
#include "cctypes.h"
#include "filefuncs.h"
#include "fslfuncs.h"
#include "newimage/newimageall.h"

extern "C" __declspec(dllexport) int _stdcall fslstats(char *CmdLn);
float mult = 1;
float thr = 0.5;

NEWIMAGE::volume<float> entrada;
read_volume(entrada, string("\\\\10.36.4.201\\dados1\\PROJETOS\\PRJ1305_UMBRELLA_TOC\\03_PROCS\\PROC_DATA\\SEARCHLIGHT\\SVC_linear_R10_CV36_PAIRED_18\\BROCCOLI\\SVM_123vx\\CULPA_NEUTRO\\CULPA_NEUTRO_4D_acc.nii"));
entrada *= mult;
NEWIMAGE::tfce(entrada, 2, 0.5, 26, thr * mult, 0);

float maxval = 0;
float maxvalReal = entrada.max();

for (int t = 0; t < 10; t++)
{
NEWIMAGE::volume<float> tempvol;
char filename[1000];

sprintf(filename, "\\\\10.36.4.201\\dados1\\PROJETOS\\PRJ1305_UMBRELLA_TOC\\03_PROCS\\PROC_DATA\\SEARCHLIGHT\\SVC_linear_R10_CV36_PAIRED_18\\BROCCOLI\\SVM_123vx\\CULPA_NEUTRO\\PERM\\CULPA_NEUTRO_4D_acc_%.5d.nii", t + 1);
read_volume(tempvol, string(filename));
tempvol *= mult;

NEWIMAGE::tfce(tempvol, 2, 0.5, 26, thr * mult, 0);

float actualMax = tempvol.max();
if (maxval < actualMax) maxval = actualMax;

printf("%d - %f\n", t + 1, actualMax);
}
printf("Perm - %f - Real - %f\n", maxval, maxvalReal);
return 0;
*/
