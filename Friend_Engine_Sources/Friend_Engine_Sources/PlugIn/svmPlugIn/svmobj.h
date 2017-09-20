#include "svm.h"
#include "masks.h"
#include "logObject.h"

// Class to relocate the mean of a projection time serie
class AdaptingSVM
{
	float positiveClass, negativeClass;
	IncrementalStats incrementalMean;

public:
	void initialize(char *model);
	void initialize(struct svm_model *model);
	void adaptResult(float &classNumber, float &projection, int adapt = 1);
};

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

	// variable the holds the logging system
	LogObject *logObject;
#ifndef __GNUC__
   SVMObj () { predict_probability=0; max_nr_attr = 64; line = NULL; };
#else //The GCC way
   SVMObj () { predict_probability=0; max_nr_attr = 64; line = NULL; };
#endif
   virtual ~SVMObj() { };
     
};

// get a svm score based on a vector sample. this functions returns the prediction class as a return value.
double svmSamplePredict(const svm_model *model, const svm_node *sample, double &score);
void initParam(struct svm_parameter &params);
void initProblem(struct svm_problem &problem, int rows, int cols, int allocateX = 1);

// deallocates problem previously initialized
void deallocateProblem(struct svm_problem &problem);


