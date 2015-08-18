#include <stdlib.h>
#include <string.h>

char* deblank(char* input)
{
	char *output = input;
	int i, j;
	for (i = 0, j = 0; i<strlen(input); i++, j++)
	{
		if (input[i] != ' ')
			output[j] = input[i];
		else
			j--;
	}
	output[j] = '\0';
	return output;
}

// parsers a string in argc and argv variables
void parser(const char *CmdLn, int &argc, char** &argv, char separator)
{
   char *auxC=NULL, *token, *unparsedString;
   unsigned int maxLen = 0, i;
   
   // just copying the CmdLn variable
   auxC = (char *) malloc(strlen(CmdLn)+5);
   strcpy(auxC, CmdLn);

   argc = 0;

   for(; auxC[0]== separator && auxC[0]!=0; auxC++); // jumping separator
   unparsedString=auxC;

   // first pass just to count tokens
   while (unparsedString != NULL)
   {
        argc++;
        if (unparsedString[0] != '"')
           token = strchr(unparsedString, separator);
        else
           token = strchr(++unparsedString, '"');

		if (token != NULL)
        {
           if (token-unparsedString > maxLen) maxLen = token-unparsedString;
           token++;
           for(; token[0]==separator && token[0]!=0; token++); //jumping separator
        }
        else
           if (strlen(unparsedString) > maxLen) maxLen = strlen(unparsedString);
        unparsedString=token;
        if (unparsedString != NULL)
           if (strlen(unparsedString) == 0) unparsedString=NULL;
   }

   // allocating memory
   argv = (char **) malloc(argc * sizeof(char *));

   unparsedString=auxC;

   // second pass to retain tokens
   i = 0;
   while (unparsedString != NULL)
   {
	    if (unparsedString[0] != '"')
           token = strchr(unparsedString, separator);
        else
           token = strchr(++unparsedString, '"');

		if (token != NULL)
        {
           token[0] = '\0';
           argv[i] = (char *) malloc(maxLen + 1);
           strcpy(argv[i], unparsedString);

           token++;
           for(; token[0]== separator && token[0]!=0; token++); // jumping separator
           i++;
        }
        else if (unparsedString != NULL)
        {
           argv[i] = (char *) malloc(maxLen + 1);
           strcpy(argv[i], unparsedString);
        }
        unparsedString=token;
        if (unparsedString != NULL)
           if (strlen(unparsedString) == 0) unparsedString=NULL;
   }
   free(auxC);
}

// frees the argv variable
void freeparser(int &argc, char** &argv)
{
  for(int i=0; i<argc;i++) free(argv[i]);
  free(argv);
}
