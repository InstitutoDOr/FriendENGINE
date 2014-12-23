#include <stdio.h>
#include "engine.h"
#include "cctypes.h"
#include <sstream>
#include <cstring>
#include "filefuncs.h"
#include <ctime>
#ifdef WINDOWS
#include <direct.h>
#endif

// Usage : engine PORT_NUMBER
// where PORT_NUMBER is the TCP/IP listening port that FRIEND Engine uses to
// receive the commands from the FRONTEND. The core functionality is in the
// FRIEND Engine object

void MniToSubject(char *betRFI, char *mniTemplate, char *mniStandard, char* output, char *prefix);

int main(int argc, char* argv[])
{
   friendEngine app;
   char port[20] = "5678";

   MniToSubject("E:\\projetos\\NFB_VR\\RFI.nii", "E:\\projetos\\NFB_VR\\Mask_VR_ROI.nii", "E:\\projetos\\NFB_VR\\MNI152_T1_1mm_brain.nii", "E:\\projetos\\NFB_VR\\rfimask.nii", NULL);
   extractFilePath(argv[0], app.workingDir);
   if (argc >1) strcpy(port, argv[1]);
   // entry function for the object. This starts the FRIEND Engine server
   app.server((BYTE *)port);
   printf("Leaving...\n");
   return 0;
}
