#include <iostream>
#include <armadillo>
#include "newimage/newimageall.h"

using namespace std;
using namespace arma;
using namespace NEWIMAGE;

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html

typedef std::vector<double> doubleVect;

void Matrix2Mat(Matrix &a, mat &b);
void Mat2Matrix(mat &a, Matrix &b);
void znormalization_arma(mat& matrix, const int dim = 0);
void ridgeArma(mat &regressors, mat &data, doubleVect &alphas, mat &betas, int normalpha = 0, double singcutoff = 1e-10);
void ridgeArma2(mat &regressors, mat &data, doubleVect &alphas, mat &betas, mat &corr, int normalpha = 0, double singcutoff = 1e-10);
void ridge_corr(mat &regressorsTrain, mat &regressorsTest, mat &dataTrain, mat &dataTest, doubleVect &alphas, mat &correlations, int normalpha = 0, double singcutoff = 1e-10);
void calculateBetas(mat &regressors, mat &data, doubleVect &alphas, mat &betas, int normalpha = 0, double singcutoff = 1e-10);
void crossValidation(mat & regressors, mat &data, doubleVect &alphas, doubleVect &bestAlphas, int nfolds=10, unsigned int seed=0);
