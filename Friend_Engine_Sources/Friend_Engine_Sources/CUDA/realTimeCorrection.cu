#include "RealTimeCorrection.cuh"
#include "legendre.h"

void printMatrix(int ncols, int n, const float*A, int lda, const char* name)
{
	for (int row = 0; row < ncols; row++){
		for (int col = 0; col < n; col++){
			double Areg = A[row + col*lda];
			//			printf("%s(%d,%d) = %f\n", name, row + 1, col + 1, Areg);
			printf("%5.5f\t", Areg);
		}
		printf("\n");
	}
}

void printMatrixT(int n, int ncols, const float*A, int lda, const char* name)
{
	for (int row = 0; row < n; row++){
		for (int col = 0; col < ncols; col++){
			double Areg = A[row*ncols + col];
			printf("%5.5f\t", Areg);
		}
		printf("\n");
	}
}

void fprintMatrix(int ncols, int n, const float*A, int lda, const char* name)
{
	FILE *f = fopen(name, "wt+");
	for (int row = 0; row < ncols; row++){
		for (int col = 0; col < n; col++){
			double Areg = A[row + col*lda];
			fprintf(f, "%5.5f,", Areg);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void fprintMatrixT(int n, int ncols, const float*A, int lda, const char* name)
{
	FILE *f = fopen(name, "wt+");
	for (int row = 0; row < n; row++){
		for (int col = 0; col < ncols; col++){
			double Areg = A[row*ncols + col];
			fprintf(f, "%5.5f,", Areg);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

__global__ void computeTValue(int Ydim, int Ndata, float *Y, float *yret, float *contrastArray)
{
	int voxel = blockIdx.x * blockDim.x + threadIdx.x;
	if (voxel < Ydim)
	{
		float oldMean = 0, variance = 0, mean = 0;
		for (int volume = 0; volume < Ndata; volume++)
		{
			float value = Y[voxel*Ndata + volume];

			oldMean = mean;
			mean = mean + (value - oldMean) / (volume + 1);

			variance += (value - oldMean)*(value - mean);
		}
		variance /= (Ndata - 1);
		yret[voxel] = yret[voxel] / sqrt(variance * contrastArray[0]);
	}

}

__global__ void zeroConstantRow(int Ydim, int line, int xdim_cur, float *XtY)
{
	int voxel = blockIdx.x * blockDim.x + threadIdx.x;
	if (voxel < Ydim)
		XtY[line + voxel * xdim_cur] = 0;
}

/***************/
/* COPY KERNEL */
/***************/
__global__ void copy_kernel(const float * __restrict d_in1, float * __restrict d_out1, const float * __restrict d_in2, float * __restrict d_out2, const int M, const int N) {

	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	const int j = blockIdx.y * blockDim.y + threadIdx.y;

	if ((i < N) && (j < N)) {
		d_out1[j * N + i] = d_in1[j * M + i];
		d_out2[j * N + i] = d_in2[j * M + i];
	}
}

/* kernel function to normalize X */
__global__ void normalizeX(int m, int n, float *dev_X)
{
	/* Normalize devX to row mean=0, SD=100
	* Input:
	* 	m,n; row and column numbers of dev_X
	* 	dev_X; float pointer to device matrix
	* Output:
	* 	dev_X is overwritten by the normalized value
	*/

	float sum, ssum; // sum and squared sum
	float mean, sd; // mean and standard deviation
	int xi; // row index of X
	int j; // column index of X
	float x;

	xi = blockIdx.x * blockDim.x + threadIdx.x;
	if (xi >= m)
		return;

	sum = 0.0;
	ssum = 0.0;
	for (j = 0; j < n; j++) {
		x = *(dev_X + xi + (j * m));
		sum += x;
		ssum += x * x;
	}

	mean = sum / n;
	sd = sqrt(ssum / n - mean * mean);

	if (sd < FLT_EPSILON)
		return;

	for (j = 0; j < n; j++) {
		// subtract mean
		*(dev_X + xi + (j * m)) -= mean;
		// divide by sd x 100
		/* Calculation error could be large for small X values.
		* For keeping accuracy, X values are scaled to SD=100
		*/
		*(dev_X + xi + (j * m)) *= (100.0 / sd);
	}
};

/* kernel function to scaling Y */
__global__ void scalingY(int m, int n, float *dev_Ymeanscal, float *dev_Y) {
	/* Scaling devX to percent change
	* Input:
	* 	m,n; row and column numbers of dev_Y
	* 	dev_Ymeanscal; float pointer to device matrix of scaling factor
	* 	dev_Y; float pointer to device matrix of scaled value
	* Output:
	* 	dev_Y is overwritten by the scaled value
	*/

	int yi; // row index of Y
	int j; // column index of Y

	yi = blockIdx.x * blockDim.x + threadIdx.x;
	if (yi >= m)
		return;

	if (*(dev_Ymeanscal + yi) == 0.0) {
		float sum = 0.0;
		float mean;
		for (j = 0; j < n; j++)
			sum += *(dev_Y + yi + (j * m));

		mean = sum / float(n);
		if (mean >= FLT_EPSILON)
			*(dev_Ymeanscal + yi) = 100.0 / mean;
	}

	for (j = 0; j < n; j++)
		*(dev_Y + yi + (j * m)) *= *(dev_Ymeanscal + yi);
};


int RealTimeCorrection::rt_glm_cuda_init(int xdim, int ydim, int nmax, int rt_reg_ricor, int poly_order, float tr_sec)
{
	/* Initialize variables and keep memory
	* Input:
	* 	xdim, ydim; data dimensionality of X (regressor) and Y (regressed data)
	* 	nmax; maximum number of samples for GLM
	* 	rt_reg_ricor; bool flag for regress RVT and ricor
	* 	poly_order; maximum order of polynomial regressor
	* 	            (This may be overwritten by the estimated value with nmax)
	* 	tr_sec; Duration of one volume (sec)
	* 	NormX; bool flag for normalizing X (mean=0, std=1)
	* 	ScaleY; bool flag for scaling Y to percent change
	*/

	cudaError_t cudaStat;
	/*--- Set parameters on global variables ---*/
	Ndata = 0;
	po_current = 0;
	Ydim = ydim;
	Nmax = nmax;
	TR_sec = tr_sec;
	ricor = rt_reg_ricor;
	data_xd = xdim;

	// Set regressor dimension
	Xdim = data_xd;

	// Add number of physiological noise regressors
	if (ricor) {
		row_ricor = Xdim;
		Xdim += 13; // Add RVT*5 + Resp*4 + ECG*4
	}

	// Update max_poly_order
	float len_sec = TR_sec * Nmax;
	int po = 1 + (int)(len_sec / 150.0); // polynomial order
	if (poly_order > po)
		max_poly_order = po;
	else
		max_poly_order = poly_order;

	// Add number of polynomial regressors to Xdim
	row_poly = Xdim;
	if (max_poly_order > 0)
		Xdim += max_poly_order + 1; // +1 is the baseline regressor (poly_order = 0)

	/*--- Allocate memory ---*/
	//Initialize pointers for conditional free
	X_buff = NULL;
	X = NULL;
	Y_buff = NULL;
	Y = NULL;
	Ymeanscal = NULL;
	XtX = NULL;
	iXtX = NULL;
	XtXa = NULL;
	iXtXa = NULL;
	XtY = NULL;
	info = NULL;
	contrast = NULL;

	cudaStat = cudaMalloc((void**)&info, sizeof(int));
	CUDA_CHECK_RETURN(cudaStat);

	// XtX
	cudaStat = cudaMalloc((void**)&XtX, Xdim * Xdim * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	cudaStat = cudaMalloc((void**)&iXtX, Xdim * Xdim * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	// XtY
	cudaStat = cudaMalloc((void**)&XtY, Xdim * Ydim * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	// XtX array
	cudaStat = cudaMalloc((void**)&XtXa, sizeof(float *));
	CUDA_CHECK_RETURN(cudaStat);

	// iXtX array
	cudaStat = cudaMalloc((void**)&iXtXa, sizeof(float *));
	CUDA_CHECK_RETURN(cudaStat);

	float **buff;
	buff = &XtX;
	cudaMemcpy(XtXa, buff, sizeof(float *), cudaMemcpyHostToDevice);

	buff = &iXtX;
	cudaMemcpy(iXtXa, buff, sizeof(float *), cudaMemcpyHostToDevice);

	// X_buff; ring buffer for X
	cudaStat = cudaMalloc((void**)&X_buff, Xdim * Nmax * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);
	// Initialize with 0
	cudaStat = cudaMemset((void*)X_buff, 0, Xdim * Nmax * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	// X
	cudaStat = cudaMalloc((void**)&X, Xdim * Nmax * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);
	// Initialize with 0
	cudaStat = cudaMemset((void*)X, 0, Xdim * Nmax * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	// Y_buff; ring buffer for Y
	cudaStat = cudaMalloc((void**)&Y_buff, Ydim * Nmax * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);
	// Initialize with 0
	cudaStat = cudaMemset((void*)Y_buff, 0, Ydim * Nmax * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	// Y
	cudaStat = cudaMalloc((void**)&Y, Ydim * Nmax * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);
	// Initialize with 0
	cudaStat = cudaMemset((void*)Y, 0, Ydim * Nmax * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	residuals = (float *) malloc(sizeof(float) * Nmax * Ydim);

	// Ymeanscal
	cudaStat = cudaMalloc((void**)&Ymeanscal, Ydim * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	// Initialize with 0
	cudaStat = cudaMemset((void*)Ymeanscal, 0, Ydim * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	// yret; for returning processed value vector of Y
	cudaStat = cudaMalloc((void**)&yret, Ydim * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	if (max_poly_order > 0)
	{
		/*-- Set baseline value in X --*/
		float one = 1.0;
		int n;
		for (n = 0; n < Nmax; n++)
		{
			cudaStat = cudaMemcpy(X_buff + row_poly + n * Xdim, &one, sizeof(float),
				cudaMemcpyHostToDevice);
			CUDA_CHECK_RETURN(cudaStat);
		}
	}
	return 0;
}

int RealTimeCorrection::rt_glm_cuda_append_data(int t, float *x, float *y) {
	/* Add data vector x, y to X, Y
	* X and Y are saved in row-order for convenience of online update
	* Data size N is updated to t+1
	* Input:
	* 	t: point to insert data (0 base);
	*  x: pointer to vector of design matrix (Xdim length)
	*  y: pointer to vector of signal matrix (Ydim length)
	*/

	cudaError_t cudaStat;

	Ndata = t + 1; //update the number of data

	/*--- Set data write position in the double buffer ---*/
	int buf_write_pos;
	buf_write_pos = t % Nmax;

	// Copy x at the buf_write_pos of dev_X_buff
	cudaStat = cudaMemcpy(X_buff + buf_write_pos * Xdim, x,
		data_xd * sizeof(float), cudaMemcpyHostToDevice);
	CUDA_CHECK_RETURN(cudaStat);

	// Copy y at the buf_end of dev_Y_buff
	cudaStat = cudaMemcpy(Y_buff + (buf_write_pos * Ydim), y,
		Ydim * sizeof(float), cudaMemcpyHostToDevice);
	CUDA_CHECK_RETURN(cudaStat);

	return 0;
}

int RealTimeCorrection::setContrastNumber(int count)
{
	cudaError_t cudaStat;
	contrastCount = count;
	cudaStat = cudaMalloc((void**)&contrast,  contrastCount * Xdim * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);

	cudaStat = cudaMemset((void*)contrast, 0, contrastCount * Xdim * sizeof(float));
	CUDA_CHECK_RETURN(cudaStat);
	return 0;
}

int RealTimeCorrection::addContrast(float *contrst, int index, int len)
{
	cudaError_t cudaStat; 
	cudaStat = cudaMemcpy(contrast + index * Xdim, contrst, len * sizeof(float), cudaMemcpyHostToDevice);
	CUDA_CHECK_RETURN(cudaStat);
	return 0;
}

int RealTimeCorrection::rt_glm_cuda_update_ricor(float *X_ricor, int len) {
	/* Update ricor regressors with len length
	* Input:
	*     X_ricor; pointer to that ricor regressor data (13 row * len column)
	*     len; length of regressor
	*/

	// If len is shorter than Ndata (if Ndata<Nmax) or Nmax, return error
	int data_len;
	if (Ndata < Nmax)
		data_len = Ndata;
	else
		data_len = Nmax;

	if (len < data_len) {
		fprintf(stderr, "RT: rt_glm_cuda, length of ricor is too short.\n");
		return -1;
	}

	// Copy ricor regressors
	int n;
	cudaError_t cudaStat;
	for (n = 0; n < data_len; n++) {
		cudaStat = cudaMemcpy(X_buff + row_ricor + n * Xdim, X_ricor + n * 13,
			sizeof(float) * 13, cudaMemcpyHostToDevice);
		CUDA_CHECK_RETURN(cudaStat);
	}

	return 0;
}

int RealTimeCorrection::rt_glm_cuda_update_polyreg()
{
	/* Update polynomial regressors up to the current saved data length (Ndata)
	*/

	// If Ndata is not enough even for linear trend, return
	if (Ndata < 1)
		return -1;

	// Set length of regressor
	int polyreg_len;
	if (Ndata < Nmax)
		polyreg_len = Ndata;
	else
		polyreg_len = Nmax;

	// Set polynomial order
	if (polyreg_len < Nmax) {
		float len_sec = TR_sec * polyreg_len;
		po_current = 1 + (int)floor(len_sec / 150.0); // polynomial order
		if (po_current > max_poly_order)
			po_current = max_poly_order;
	}
	else
		po_current = max_poly_order;

	// Denominator for mapping 0..N_polyreg-1 into -1..1
	double aa = 2.0 / (polyreg_len - 1);
	double xx;

	// Append polynomial regressors
	int po, n;
	float poly_x;
	cudaError_t cudaStat;
	for (po = 1; po < max_poly_order + 1; po++) {
		for (n = 0; n < polyreg_len; n++) {
			if (po <= po_current) {
				xx = aa * n - 1.0;
				poly_x = (float)legendre(xx, po);
			}
			else
				poly_x = 0.0;

			cudaStat = cudaMemcpy(X_buff + row_poly + po + n * Xdim, &poly_x,
				sizeof(float), cudaMemcpyHostToDevice);
			CUDA_CHECK_RETURN(cudaStat);
		}
	}

	return 0;
}

int RealTimeCorrection::rt_glm_cuda_get_residual(int t, int normX, int scaleY, float *YRet)
{
	/* Get residual vector of Y at t after regressing out X from Y
	* Input:
	* 	t: point of Y to get result (0 base)
	* 	normX; bool flag for normalizing X
	* 	scaleY; bool flag for scaling Y to percent change
	* Output:
	* 	yret: pointer to return vector (sizeof(float)*Ydim memory must be
	* 	      allocated)
	*/

	cudaError_t cudaStat;
	char filename[100];
	int nthreads = 512;

	/*-- Check validity of t --*/
	if (t >= Ndata || t < Ndata - Nmax) { // Error in t
		fprintf(stderr, "RT: Error (rt_glm_cuda_get_residual); ");
		fprintf(stderr, "Time point %d is not in buffer ", t);
		if (Ndata < Nmax)
			fprintf(stderr, "(0-%d)\n", Ndata - 1);
		else
			fprintf(stderr, "(%d-%d)\n", Ndata - Nmax, Ndata - 1);

		memset(YRet, 0, sizeof(float) * Ydim);
		return -1;
	}

	/*-- Set length of data for regression --*/
	int N_reg; // Length of data for regression
	if (Ndata < Nmax)
		N_reg = Ndata;
	else
		N_reg = Nmax;

	// Current Xdim (exclude unset polynomial)
	int xdim_cur = Xdim - max_poly_order + po_current;

	cublasStatus_t cublas_status;

	cublas_status = cublasSgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, Ndata, Xdim, &one, X_buff, Xdim, &zero, X_buff, Ndata, X, Ndata);
	cublas_status = cublasSgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, Ndata, Ydim, &one, Y_buff, Ydim, &zero, Y_buff, Ndata, Y, Ndata);

	if (1) //(outputFiles)
	{
		float *A = (float*)malloc(sizeof(float) * Xdim * Ndata);
		cudaMemcpy(A, X, sizeof(float) * Ndata * Xdim, cudaMemcpyDeviceToHost);

		sprintf(filename, "E:\\debug\\A_%d.csv", Ndata);
		fprintMatrix(Ndata, xdim_cur, A, Ndata, filename);

		free(A);
	}

	if (outputFiles)
	{
		float *A = (float*)malloc(sizeof(float) * Xdim * Ndata);
		float *B = (float*)malloc(sizeof(float) * Ydim * Ndata);
		cudaMemcpy(A, X, sizeof(float) * Ndata * Xdim, cudaMemcpyDeviceToHost);
		cudaMemcpy(B, Y, sizeof(float) * Ndata * Ydim, cudaMemcpyDeviceToHost);

		sprintf(filename, "E:\\debug\\A_%d.csv", Ndata);
		fprintMatrix(Ndata, Xdim, A, Ndata, filename);

		sprintf(filename, "E:\\debug\\B_%d.csv", Ndata);
		fprintMatrix(Ndata, Ydim, B, Ndata, filename);
		free(A);
		free(B);
	}

	dim3 threads = dim3(CUDANTHREAD, 1);
	dim3 blocks;

	/*-- Normalize X --*/
	if (normX == 1) {
		blocks = dim3(Xdim / threads.x, 1);
		normalizeX << <blocks, threads >> > (Xdim, N_reg, X);
	}

	/*-- Scale Y --*/
	if (scaleY == 1) {
		blocks = dim3(Ydim / threads.x, 1);
		scalingY << <blocks, threads >> >(Ydim, N_reg, Ymeanscal, Y);
	}

	cublas_status = cublasSgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, xdim_cur, xdim_cur, Ndata, &one, X, Ndata, X, Ndata, &zero, XtX, xdim_cur);
	assert(CUBLAS_STATUS_SUCCESS == cublas_status);

	cublas_status = cublasSgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, xdim_cur, Ydim, Ndata, &one, X, Ndata, Y, Ndata, &zero, XtY, xdim_cur);
	assert(CUBLAS_STATUS_SUCCESS == cublas_status);

	if (outputFiles)
	{
		float *A = (float*)malloc(sizeof(float) * Xdim * Xdim);
		float *B = (float*)malloc(sizeof(float) * Ydim * Xdim);
		cudaMemcpy(A, XtX, sizeof(float) * xdim_cur * xdim_cur, cudaMemcpyDeviceToHost);
		cudaMemcpy(B, XtY, sizeof(float) * xdim_cur * Ydim, cudaMemcpyDeviceToHost);

		sprintf(filename, "E:\\debug\\XtX_%d.csv", Ndata);
		fprintMatrix(xdim_cur, xdim_cur, A, xdim_cur, filename);

		sprintf(filename, "E:\\debug\\XtY_%d.csv", Ndata);
		fprintMatrix(xdim_cur, Ydim, B, xdim_cur, filename);
		free(A);
		free(B);
	}

	cublas_status = cublasSmatinvBatched(handle, xdim_cur, (const float **)XtXa, xdim_cur, iXtXa, xdim_cur, info, 1);
	assert(CUBLAS_STATUS_SUCCESS == cublas_status);

	if (outputFiles)
	{
		float *A = (float*)malloc(sizeof(float) * Xdim * Xdim);
		cudaMemcpy(A, iXtX, sizeof(float) * xdim_cur * xdim_cur, cudaMemcpyDeviceToHost);

		sprintf(filename, "E:\\debug\\iXtX_%d.csv", Ndata);
		fprintMatrix(xdim_cur, xdim_cur, A, xdim_cur, filename);

		free(A);
	}

	cublas_status = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, xdim_cur, Ydim, xdim_cur, &one, iXtX, xdim_cur, XtY, xdim_cur, &zero, XtY, xdim_cur);
	assert(CUBLAS_STATUS_SUCCESS == cublas_status);

	if (outputFiles)
	{
		float *XCalc = (float*)malloc(sizeof(float) * Xdim * Ydim);
		cudaStat = cudaMemcpy(XCalc, XtY, Ydim * Xdim * sizeof(float), cudaMemcpyDeviceToHost);
		sprintf(filename, "E:\\debug\\X_%d.csv", Ndata);
		fprintMatrix(xdim_cur, Ydim, XCalc, xdim_cur, filename);
		free(XCalc);
	}

	if (calculateStatistics)
	{
		// residuals in Y
		cublas_status = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, Ndata, Ydim, xdim_cur, &alpha, X, Ndata, XtY, xdim_cur, &beta, Y, Ndata);
		assert(CUBLAS_STATUS_SUCCESS == cublas_status);

		// cT * Betas
		cublas_status = cublasSgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, contrastCount, Ydim, xdim_cur, &one, contrast, xdim_cur, XtY, xdim_cur, &zero, yret, contrastCount);
		assert(CUBLAS_STATUS_SUCCESS == cublas_status);

		if (outputFiles)
		{
			float *A = (float*)malloc(sizeof(float) * contrastCount * Ydim);
			float *B = (float*)malloc(sizeof(float) * Ydim * Ndata);
			cudaMemcpy(A, yret, sizeof(float) * contrastCount * Ydim, cudaMemcpyDeviceToHost);
			cudaMemcpy(B, Y, sizeof(float) * Ndata * Ydim, cudaMemcpyDeviceToHost);

			sprintf(filename, "E:\\debug\\CTBetas_%d.csv", Ndata);
			fprintMatrix(contrastCount, Ydim, A, contrastCount, filename);

			sprintf(filename, "E:\\debug\\Residuals_%d.csv", Ndata);
			fprintMatrix(Ndata, Ydim, B, Ndata, filename);
			free(A);
			free(B);
		}

		// cT * iXtX
		cublas_status = cublasSgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, contrastCount, xdim_cur, xdim_cur, &one, contrast, xdim_cur, iXtX, xdim_cur, &zero, XtX, xdim_cur);
		assert(CUBLAS_STATUS_SUCCESS == cublas_status);

		// anterior * c
		cublas_status = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, contrastCount, xdim_cur, xdim_cur, &one, XtX, xdim_cur, contrast, xdim_cur, &zero, XtX, contrastCount);
		assert(CUBLAS_STATUS_SUCCESS == cublas_status);

		computeTValue <<<(Ydim + nthreads-1) / nthreads, nthreads>>> (Ydim, Ndata, Y, yret, XtX);

		cudaStat = cudaMemcpy(YRet, yret, Ydim * sizeof(float), cudaMemcpyDeviceToHost);
		CUDA_CHECK_RETURN(cudaStat);
	}
	else
	{
		// yret = B(t) - A(t)*X
		cudaMemcpy(yret, Y_buff + ((Ndata - 1) * Ydim), sizeof(float) * Ydim, cudaMemcpyDeviceToDevice);

		zeroConstantRow << <(Ydim + nthreads - 1) / nthreads, nthreads >> > (Ydim, row_poly, xdim_cur, XtY);
		if (outputFiles)
		{
			float *XCalc = (float*)malloc(sizeof(float) * Xdim * Ydim);
			cudaStat = cudaMemcpy(XCalc, XtY, Ydim * Xdim * sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(filename, "E:\\debug\\Xdp_%d.csv", Ndata);
			fprintMatrix(xdim_cur, Ydim, XCalc, xdim_cur, filename);
			free(XCalc);
		}

		cublas_status = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 1, Ydim, xdim_cur, &alpha, X_buff + ((Ndata - 1) * Xdim), 1, XtY, xdim_cur, &beta, yret, 1);
		assert(CUBLAS_STATUS_SUCCESS == cublas_status);

		// Copy back dev_yret to host Ys
		cudaStat = cudaMemcpy(YRet, yret, Ydim * sizeof(float), cudaMemcpyDeviceToHost);
		CUDA_CHECK_RETURN(cudaStat);
		if (outputFiles)
		{
			sprintf(filename, "E:\\debug\\BRTC_%d.csv", Ndata);
			fprintMatrix(1, Ydim, YRet, 1, filename);
		}

	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int RealTimeCorrection::rt_glm_cuda_free_GPU_memory()
{
	/* Free GPU memory */

	cudaError_t cudaStat;

	if (X_buff) {
		cudaStat = cudaFree(X_buff);
		CUDA_CHECK_RETURN(cudaStat);
	}

	if (Y_buff) {
		cudaStat = cudaFree(Y_buff);
		CUDA_CHECK_RETURN(cudaStat);
	}

	if (X) {
		cudaStat = cudaFree(X);
		CUDA_CHECK_RETURN(cudaStat);
	}

	if (Y) {
		cudaStat = cudaFree(Y);
		CUDA_CHECK_RETURN(cudaStat);
	}

	if (Ymeanscal) {
		cudaStat = cudaFree(Ymeanscal);
		CUDA_CHECK_RETURN(cudaStat);
	}

	if (yret) {
		cudaStat = cudaFree(yret);
		CUDA_CHECK_RETURN(cudaStat);
	}

	if (XtX) cudaFree(XtX);
	if (XtY) cudaFree(XtY);

	if (iXtX) cudaFree(iXtX);
	if (iXtXa) cudaFree(iXtXa);
	if (XtXa) cudaFree(XtXa);
	if (info) cudaFree(info);
	if (contrast) cudaFree(contrast);
	if (residuals) free(residuals);

	return 0;
}


/*----------------------------------------------------------------------------*/
int RealTimeCorrection::rt_glm_cuda_get_Ndata() {
	/* Return number of data points saved in the library */
	return Ndata;
}

/*----------------------------------------------------------------------------*/
/* Print matrix values on device (for debug) */
int RealTimeCorrection::_print_device_Mtx(int n, int m, float* dev_Mtx)
{
	float *host_Mtx = (float *)malloc(n * m * sizeof(float));

	cudaError_t cudaStat;
	cudaStat = cudaMemcpy(host_Mtx, dev_Mtx, n * m * sizeof(float),
		cudaMemcpyDeviceToHost);
	CUDA_CHECK_RETURN(cudaStat);

	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++)
			printf("%f ", host_Mtx[i + j * n]);

		printf("\n");
	}

	free(host_Mtx);
	return 0;
}

/*----------------------------------------------------------------------------*/
void RealTimeCorrection::print_X_on_GPU()
{
	int N_print;
	if (Ndata < Nmax * 2)
		N_print = Ndata;
	else
		N_print = Nmax * 2;

	_print_device_Mtx(Xdim, N_print, X);
}

void RealTimeCorrection::print_Y_on_GPU(void)
{
	int N_print;
	if (Ndata < Nmax * 2)
		N_print = Ndata;
	else
		N_print = Nmax * 2;

	_print_device_Mtx(Ydim, N_print, Y);
}

/*----------------------------------------------------------------------------*/
/* Get matrix values on device */
int RealTimeCorrection::_get_device_Mtx(int n, int m, float* dev_Mtx, float* host_Mtx)
{
	cudaError_t cudaStat;
	cudaStat = cudaMemcpy(host_Mtx, dev_Mtx, n * m * sizeof(float),
		cudaMemcpyDeviceToHost);
	CUDA_CHECK_RETURN(cudaStat);

	return 0;
}

/*----------------------------------------------------------------------------*/
int RealTimeCorrection::get_Xdim_on_GPU()
{
	return Xdim;
}

int RealTimeCorrection::get_Ydim_on_GPU()
{
	return Ydim;
}

void RealTimeCorrection::get_X_on_GPU(float* host_Mtx)
{
	_get_device_Mtx(Xdim, Nmax * 2, X, host_Mtx);
}

void RealTimeCorrection::get_Y_on_GPU(float* host_Mtx)
{
	_get_device_Mtx(Ydim, Nmax * 2, Y, host_Mtx);
}

void RealTimeCorrection::get_beta_on_GPU(float* host_Mtx)
{
	//	_get_device_Mtx(Xdim, Ydim, XtY, host_Mtx);
}

int RealTimeCorrection::init_cuda()
{
	status = cublasCreate(&handle);

	if (status != CUBLAS_STATUS_SUCCESS)
	{
		fprintf(stderr, "!!!! CUBLAS initialization error\n");
		return EXIT_FAILURE;
	}
	return 0;
}

int RealTimeCorrection::end_cuda()
{
	/* Shutdown */
	status = cublasDestroy(handle);

	// cudaDeviceReset causes the driver to clean up all state. While
	// not mandatory in normal operation, it is good practice.  It is also
	// needed to ensure correct operation when the application is being
	// profiled. Calling cudaDeviceReset causes all profile data to be
	// flushed before the application exits
	cudaDeviceReset();
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		fprintf(stderr, "!!!! shutdown error (A)\n");
		return EXIT_FAILURE;
	}
	return 0;
}

int RealTimeCorrection::canGetResidual()
{
	//int xdim_cur = Xdim - max_poly_order + po_current;
	//if (Ndata < Xdim) return 0;
	if (Ndata < 50) return 0;
	return 1;
}