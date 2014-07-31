// parsers a string CmdLn in argc and argv variables. This function is crucial in transforming programs in functions and you need to put the name of program as like in bash to the transformation works with minor changes
void parser(const char *CmdLn, int &argc, char** &argv, char separator = ' ');

// frees the argv variable
void freeparser(int &argc, char** &argv);
