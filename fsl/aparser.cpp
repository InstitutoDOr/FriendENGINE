#include <stdlib.h>
#include <string.h>


void parser(char *CmdLn, int &argc, char** &argv)
{
   char *AuxC=NULL, *Token, *Inicio;
   int MaxLen = 0, i;
   AuxC = (char *) malloc(strlen(CmdLn)+5);
   strcpy(AuxC, CmdLn);

   argc = 0;

   for(; AuxC[0]== ' ' && AuxC[0]!=NULL; AuxC++); // Pula os Espaços
   Inicio=AuxC;

   while (Inicio != NULL)
   {
        argc++;
/*		Aspa = strchr(Inicio, '"');
		Espaco = strchr(Inicio, ' ');
        if ((Aspa != NULL) && (Espaco != NULL))
		{
			if ((Aspa-Inicio) < (Espaco-Inicio)) Token = strchr(++Aspa, '"');
			else Token = Espaco;
		}
		else if (Aspa != NULL) Token = strchr(++Aspa, '"');
			 else Token = Espaco;
*/
        if (Inicio[0] != '"')
           Token = strchr(Inicio, ' ');
        else
           Token = strchr(++Inicio, '"');

		if (Token != NULL)
        {
           if (Token-Inicio > MaxLen) MaxLen = Token-Inicio;
           Token++;
           for(; Token[0]==' ' && Token[0]!=NULL; Token++); // Pula os Espaços
        }
        else
           if (strlen(Inicio) > MaxLen) MaxLen = strlen(Inicio);
        Inicio=Token;
        if (Inicio != NULL)
           if (strlen(Inicio) == 0) Inicio=NULL;
   }

   argv = (char **) malloc(argc * sizeof(char *));

   Inicio=AuxC;

   i = 0;
   while (Inicio != NULL)
   {
/*		Aspa = strchr(Inicio, '"');
		Espaco = strchr(Inicio, ' ');
        if ((Aspa != NULL) && (Espaco != NULL))
		{
			if ((Aspa-Inicio) < (Espaco-Inicio)) Token = strchr(++Aspa, '"');
			else Token = Espaco;
		}
		else if (Aspa != NULL) Token = strchr(++Aspa, '"');
			 else Token = Espaco;
*/
	    if (Inicio[0] != '"')
           Token = strchr(Inicio, ' ');
        else
           Token = strchr(++Inicio, '"');

		if (Token != NULL)
        {
           Token[0] = '\0';
           argv[i] = (char *) malloc(MaxLen + 1);
           strcpy(argv[i], Inicio);

           Token++;
           for(; Token[0]==' ' && Token[0]!=NULL; Token++); // Pula os Espaços
           i++;
        }
        else if (Inicio != NULL)
        {
           argv[i] = (char *) malloc(MaxLen + 1);
           strcpy(argv[i], Inicio);
        }
        Inicio=Token;
        if (Inicio != NULL)
           if (strlen(Inicio) == 0) Inicio=NULL;
   }
   free(AuxC);
}

void freeparser(int &argc, char** &argv)
{
  for(int i=0; i<argc;i++) free(argv[i]);
  free(argv);
}