#include "svm.h"

// class responsible for handling the libsvm functions
class SVMObj
{
    char *line;
    int max_line_len;
    struct svm_parameter param;		// set by parse_command_line
    struct svm_problem prob;		// set by read_problem
    struct svm_model *model;
    struct svm_node *x_space;
    int cross_validation;
    int nr_fold;
   
    struct svm_node *x;
    int max_nr_attr;
    int predict_probability;


    // parsers the command line string of libsvm commands
    int parse_command_line(int argc, char **argv, char *input_file_name, char *model_file_name);

    // reads a svm examples file
    int read_problem(const char *filename);

    // do the cross validation in svm_train command
    void do_cross_validation();
   
    // reads a line of a svm samples file
    char* readline(FILE *input);
   
    // internal function of predict
    void predict(FILE *input, FILE *output);

    public:
   
    // svmtrain transformed in a function
    int train(const char *cmd);
   
    // svmpredict transformed in a function
    int predict(const char *cmd);
#ifndef __GNUC__
   SVMObj () { predict_probability=0; max_nr_attr = 64; line = NULL; };
#else //The GCC way
   SVMObj () { predict_probability=0; max_nr_attr = 64; line = NULL; };
#endif
   virtual ~SVMObj() { };
     
};

