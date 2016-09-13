#include <stdio.h>
#include "engine.h"
#include "cctypes.h"
#include <sstream>
#include <cstring>
#include "fslfuncs.h"
#include "filefuncs.h"
#include "masks.h"
#include <ctime>
#ifdef WINDOWS
#include <direct.h>
#endif

   /*
   char outputFile[500], prefix[500] = "_RFI2";

   sprintf(outputFile, "E:\\projetos\\NFB_VR\\saida_RFI2.nii");
   MniToSubject("E:\\projetos\\NFB_VR\\SUBJ010\\input\\RFI_sks.nii", "E:\\projetos\\NFB_VR\\Mask_VR_ROI.nii", "E:\\projetos\\NFB_VR\\MNI152_T1_1mm_brain.nii", outputFile, prefix);
   sprintf(outputFile, "E:\\projetos\\NFB_VR\\saida_RFI3.nii");
   MniToSubject("E:\\projetos\\NFB_VR\\SUBJ010\\input\\RFI_sks.nii", "E:\\projetos\\NFB_VR\\Mask_VR_ROI_std.nii", "E:\\projetos\\NFB_VR\\MNI152_T1_1mm_brain.nii", outputFile, prefix);

   fslstats("fslstats E:\\projetos\\NFB_VR\\RFI_sks.nii -k E:\\projetos\\NFB_VR\\saida_RFI2.nii -m");
   fslstats("fslstats E:\\projetos\\NFB_VR\\RFI_sks.nii -k E:\\projetos\\NFB_VR\\saida_RFI3.nii -m");

   return 0;

   RoiMeanCalculation meanCalculation;
   meanCalculation.loadReference("E:\\projetos\\NFB_VR\\Mask_VR_ROI_std.nii");
   meanCalculation.calculateMeans("E:\\projetos\\NFB_VR\\MNI152_T1_1mm_brain.nii");
   fprintf(stderr, "Roi 1 = %f\n", meanCalculation.roiMean(0));
   fprintf(stderr, "Roi 2 = %f\n", meanCalculation.roiMean(1));

   RoiMeanCalculation meanCalculation;
   int t = 1;
   meanCalculation.loadReference("E:\\projetos\\NFB_VR\\saida_RFI2.nii");
   meanCalculation.calculateMeans("E:\\projetos\\NFB_VR\\RFI_sks.nii");
   fprintf(stderr, "%d - %f\t%f\n", t, meanCalculation.roiMean(0), meanCalculation.roiMean(1));

   meanCalculation.loadReference("E:\\projetos\\NFB_VR\\saida_RFI3.nii");
   meanCalculation.calculateMeans("E:\\projetos\\NFB_VR\\RFI_sks.nii");
   fprintf(stderr, "%d - %f\t%f\n", t, meanCalculation.roiMean(0), meanCalculation.roiMean(1));

   return 0;

//   RoiMeanCalculation meanCalculation;
   meanCalculation.loadReference("E:\\projetos\\NFB_VR\\saida_RFI2.nii");
   for (int t = 1; t <= 304; t++)
   {
	   char fname[500];
	   sprintf(fname, "E:\\projetos\\NFB_VR\\SUBJ010\\preproc_RUN01\\DRIN-%.5d_mc_g.nii", t);
	   meanCalculation.calculateMeans(fname);
	   fprintf(stderr, "%d - %f\t%f\n", t, meanCalculation.roiMean(0), meanCalculation.roiMean(1));
   }

   return 0;
   */
   /*
   if (0)
   {
	   sprintf(outputFile, "E:\\projetos\\NFB_VR\\saida_RFI2.nii");
	   MniToSubject("E:\\projetos\\NFB_VR\\SUBJ010\\input\\RFI_sks.nii", "E:\\projetos\\NFB_VR\\Mask_VR_ROI.nii", "E:\\projetos\\NFB_VR\\MNI152_T1_1mm_brain.nii", outputFile, prefix);
   }
   else
   {
	   sprintf(outputFile, "E:\\projetos\\NFB_VR\\saida_RFI2.nii");
	   MniToSubject("E:\\projetos\\NFB_VR\\SUBJ010\\input\\RFI_sks.nii", "E:\\projetos\\NFB_VR\\Mask_VR_ROI.nii", "E:\\projetos\\NFB_VR\\MNI152_T1_1mm_brain.nii", outputFile, prefix);
	   sprintf(outputFile, "E:\\projetos\\NFB_VR\\saida_RFI3.nii");
	   MniToSubject("E:\\projetos\\NFB_VR\\SUBJ010\\input\\RFI_sks.nii", "E:\\projetos\\NFB_VR\\Mask_VR_ROI_std.nii", "E:\\projetos\\NFB_VR\\MNI152_T1_1mm_brain.nii", outputFile, prefix);
   }

   return 0;
 */
