#ifndef ____cudaRealTimeCorrection__
#define ____cudaRealTimeCorrection__

/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* Includes, cuda */

#include <math.h>
#include <float.h>
#include "cuda.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "device_launch_parameters.h"

class RealTimeCorrection
{
	cublasStatus_t status;
	cublasHandle_t handle;

	int data_xd; // length of input regressor data
	int Xdim; // Dimensionality of regressor (not include polynomial)
	int Ydim; // Dimensionality of dependent variable (MRI volume)
	int Nmax; // Maximum size of samples (Nmax*2 memory space is allocated for buffering)
	int Ndata; // Current data time
	int row_ricor; // start row of ricor regressors
	int row_poly; // start row of polynomial regressors
	float alpha = -1.0f;
	float beta = 1.0f;
	float one = 1;
	float zero = 0;
	int *info;

	// RVT,RICOR
	int ricor; // bool flag for regress rvt and ricor

	// Polynomial regressor
	/* if max_poly_order==0, only baseline regressor is included */
	int max_poly_order = 0; // Maximum order of polynomial regressor
	int po_current; // Polynomial order at current data length
	float TR_sec; // Duration of one volume (sec). Used for get poly_order

	// Pointer to memory
	float *X_buff, *Y_buff; // Ring buffer for X and Y
	float *X, *Y, *Ymeanscal, *yret, *XtX, *iXtX, *XtY, **XtXa, **iXtXa;
	float *residuals;
	float *contrast;
	int contrastCount;
	float *variances;


public:
	int outputFiles = 0;
	int calculateStatistics = 1;

	int init_cuda();
	int end_cuda();
	int canGetResidual();
	int rt_glm_cuda_init(int xdim, int ydim, int nmax, int rt_reg_ricor, int poly_order, float tr_sec);
	int rt_glm_cuda_append_data(int t, float *x, float *y);
	int rt_glm_cuda_update_ricor(float *X_ricor, int len);
	int rt_glm_cuda_update_polyreg();
	int rt_glm_cuda_get_residual(int t, int normX, int scaleY, float *yret);
	int rt_glm_cuda_free_GPU_memory();
	int rt_glm_cuda_get_Ndata();
	void get_beta_on_GPU(float* host_Mtx);
	void get_Y_on_GPU(float* host_Mtx);
	void get_X_on_GPU(float* host_Mtx);
	int get_Ydim_on_GPU();
	int get_Xdim_on_GPU();
	int _get_device_Mtx(int n, int m, float* dev_Mtx, float* host_Mtx);
	void print_Y_on_GPU(void);
	void print_X_on_GPU();
	int _print_device_Mtx(int n, int m, float* dev_Mtx);
	int addContrast(float *contrst, int index, int len);
	int setContrastNumber(int count);

};

/* Utility macro function to print error and return */
#define PERROR {\
	fprintf(stderr, "Error at line %d in file %s\n",  __LINE__, __FILE__);   \
    status = cublasDestroy(handle);                                          \
	cudaDeviceReset();                                                       \
	return 1;                                                                \
};

/* Utility macro function for checking cudaError_t value */
#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		return 1;															\
					} };

#define BLOCK_SIZE 32
#ifndef CUDANTHREAD
#define CUDANTHREAD 512
#endif

void printMatrix(int ncols, int n, const float*A, int lda, const char* name);
void printMatrixT(int n, int ncols, const float*A, int lda, const char* name);
void fprintMatrix(int ncols, int n, const float*A, int lda, const char* name);
void fprintMatrixT(int n, int ncols, const float*A, int lda, const char* name);

#endif /* defined(____cudaRealTimeCorrection__) */