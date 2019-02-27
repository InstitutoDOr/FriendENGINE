#include "sgolay.h"
#include <iostream>
#include <fstream>

using namespace std;

void readVector(char *filename, float_vect &vetor);
void saveVector(char *filename, float_vect &vetor);
float_vect detrendVector(float_vect &unfiltered, int w, int degree);
