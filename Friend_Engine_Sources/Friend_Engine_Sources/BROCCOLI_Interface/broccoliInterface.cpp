#include "broccoliInterface.h"

int             IMAGE_REGISTRATION_FILTER_SIZE = 7;
int             MOTION_CORRECTION_FILTER_SIZE = 7;
int             NUMBER_OF_ITERATIONS_FOR_MOTION_CORRECTION = 5;
int             NUMBER_OF_MOTION_CORRECTION_PARAMETERS = 6;
int             NUMBER_OF_AFFINE_IMAGE_REGISTRATION_PARAMETERS = 12;
int             NUMBER_OF_ITERATIONS_FOR_LINEAR_IMAGE_REGISTRATION = 10;
int             NUMBER_OF_ITERATIONS_FOR_NONLINEAR_IMAGE_REGISTRATION = 10;
int             COARSEST_SCALE = 4;
int             MM_T1_Z_CUT = 0;
bool            DEBUG = false;
const char*     FILENAME_EXTENSION = "_MNI";
bool            PRINT = true;
bool            WRITE_TRANSFORMATION_MATRIX = false;
bool            WRITE_DISPLACEMENT_FIELD = false;
bool			WRITE_INTERPOLATED = false;
bool			CHANGE_OUTPUT_FILENAME = false;
float			SIGMA = 5.0f;
bool			MASK = false;
bool			MASK_ORIGINAL = false;
int             OPENCL_PLATFORM = 0;
int             OPENCL_DEVICE = 0;
bool			VERBOS = false;

size_t IMAGE_REGISTRATION_PARAMETERS_SIZE = NUMBER_OF_AFFINE_IMAGE_REGISTRATION_PARAMETERS * sizeof(float);
size_t FILTER_SIZE = IMAGE_REGISTRATION_FILTER_SIZE * IMAGE_REGISTRATION_FILTER_SIZE * IMAGE_REGISTRATION_FILTER_SIZE * sizeof(float);
size_t FILTER_SIZE2 = IMAGE_REGISTRATION_FILTER_SIZE * IMAGE_REGISTRATION_FILTER_SIZE * IMAGE_REGISTRATION_FILTER_SIZE * sizeof(cl_float2);
size_t PROJECTION_TENSOR_SIZE = NUMBER_OF_FILTERS_FOR_NONLINEAR_REGISTRATION * sizeof(float);
size_t FILTER_DIRECTIONS_SIZE = NUMBER_OF_FILTERS_FOR_NONLINEAR_REGISTRATION * sizeof(float);

#ifdef WINDOWS
int getMilliCount()
{
	timeb tb;
	ftime(&tb);
	int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
	return nCount;
}

int getMilliSpan(int nTimeStart)
{
	int nSpan = getMilliCount() - nTimeStart;
	if (nSpan < 0)
		nSpan += 0x100000 * 1000;
	return nSpan;
}
#endif

double GetWallTime()
{

#ifdef WINDOWS
	return (double)getMilliCount() / 1000.0;
#else
	struct timeval time;
	if (gettimeofday(&time, NULL))
	{
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
#endif
}

void ReadBinaryFile(float* pointer, int size, const char* filename)
{
	FILE *fp = NULL;
	fp = fopen(filename, "rb");

	if (fp != NULL)
	{
		fread(pointer, sizeof(float), size, fp);
		fclose(fp);
	}
}

void vector2volume(volume4D<float> &volume, float *volumeVector)
{
	int64_t sizeVolume = volume.xsize() * volume.ysize() * volume.zsize();

	int voxel = 0;
	for (int t = 0; t < volume.tsize(); t++)
	{
		float *volumeData = (float *)volume[t].fbegin();
		for (int j = 0; j < sizeVolume; j++)
			volumeData[j] = volumeVector[voxel++];
	}
}

void vector2volume(volume<float> &volume, float *volumeVector)
{
	int64_t sizeVolume = volume.xsize() * volume.ysize() * volume.zsize();

	int voxel = 0;
    float *volumeData = (float *)volume.fbegin();
	for (int j = 0; j < sizeVolume; j++)
		volumeData[j] = volumeVector[voxel++];
}

int64_t volume2vector(volume4D<float> &volume, float *&volumeVector)
{
	int64_t sizeVolume = volume.xsize() * volume.ysize() * volume.zsize();

	volumeVector = (float *) malloc(sizeVolume * sizeof(float));

	float maxV = volume.max();
	int voxel = 0;
	for (int t = 0; t < volume.tsize(); t++)
	{
		float *volumeData = (float *) volume[t].fbegin();
		for (int j = 0; j < sizeVolume; j++)
			volumeVector[voxel++] = volumeData[j];
	}

	return sizeVolume * sizeof(float);
}

void copyVolume2vector(volume<float> &volume, float *&volumeVector)
{
	int64_t sizeVolume = volume.xsize() * volume.ysize() * volume.zsize();

	float *volumeData = (float *)volume.fbegin();
	for (int t = 0; t < sizeVolume; t++)
		volumeVector[t] = volumeData[t];
}

int64_t volume2vector(volume<float> &volume, float *&volumeVector)
{
	int64_t sizeVolume = volume.xsize() * volume.ysize() * volume.zsize();
	volumeVector = (float *)malloc(sizeVolume * sizeof(float));

	copyVolume2vector(volume, volumeVector);

	return sizeVolume * sizeof(float);
}

void BROCCOLILinearRegistration(char *T1file, char *MNIfile, char *outputFile)
{
	BROCCOLI_LIB BROCCOLI(OPENCL_PLATFORM, OPENCL_DEVICE, 2, VERBOS); // 2 = Bash wrapper
	if (BROCCOLI.GetOpenCLInitiated())
	{
		fprintf(stderr, "BROCCOLI correctly initialized !!!!\n");

		volume4D<float> T1volume;
		fprintf(stderr, "Loading T1 volume\n");
		load_volume4D(T1volume, string(T1file));

		volume<float> mniVolume;
		fprintf(stderr, "Loading MNI volume\n");
		load_volume(mniVolume, string(MNIfile));

		volume4D<float> output(mniVolume.xsize(), mniVolume.ysize(), mniVolume.zsize(), T1volume.tsize());
		for (int t = 0; t < T1volume.tsize(); t++)
			output[t].copyproperties(mniVolume);

		fprintf(stderr, "Preparing output\n");

		// Run the actual registration
		fprintf(stderr, "Running registration\n");
		prepareBROCCOLIRegistration(BROCCOLI, T1volume, mniVolume, output);
		save_volume4D(output, string(outputFile));
	}
}

void prepareBROCCOLISmoothing(BROCCOLI_LIB &BROCCOLI, volume<float> &inputData, float fwhm, float *&fMRI_Volume, float *&h_Certainty)
{
	// Set all necessary pointers and values
	int64_t allocatedHostMemory = 0;

	allocatedHostMemory += volume2vector(inputData, fMRI_Volume);
	allocatedHostMemory += (1 + 2) * volume2vector(inputData, h_Certainty);

	BROCCOLI.SetInputfMRIVolumes(fMRI_Volume);
	BROCCOLI.SetAutoMask(false);
	BROCCOLI.SetInputCertainty(h_Certainty);

	BROCCOLI.SetEPISmoothingAmount(fwhm);
	BROCCOLI.SetAllocatedHostMemory(allocatedHostMemory);

	BROCCOLI.SetEPIWidth(inputData.xsize());
	BROCCOLI.SetEPIHeight(inputData.ysize());
	BROCCOLI.SetEPIDepth(inputData.zsize());
	BROCCOLI.SetEPITimepoints(1);

	BROCCOLI.SetEPIVoxelSizeX(inputData.xdim());
	BROCCOLI.SetEPIVoxelSizeY(inputData.ydim());
	BROCCOLI.SetEPIVoxelSizeZ(inputData.zdim());

	// Run the actual slice timing correction
	BROCCOLI.PerformSmoothingNormalizedHostWrapper();
}

void prepareBROCCOLIRegistration(BROCCOLI_LIB &BROCCOLI, volume4D<float> &T1, volume<float> &mniVolume, volume4D<float> &output)
{
	float *h_T1_Volume, *h_MNI_Volume;
	float *h_Interpolated_T1_Volume, *h_Aligned_T1_Volume;
	broccoliLinearRegistrationFilters filter;
	float h_Registration_Parameters[12];
	int64_t allocatedHostMemory = 0;

	allocatedHostMemory += volume2vector(T1, h_T1_Volume);
	allocatedHostMemory += (1 + 2) * volume2vector(mniVolume, h_MNI_Volume);

	// Set all necessary pointers and values
	BROCCOLI.SetInputT1Volume(h_T1_Volume);
	BROCCOLI.SetInputMNIBrainVolume(h_MNI_Volume);

	BROCCOLI.SetAllocatedHostMemory(allocatedHostMemory);

	BROCCOLI.SetT1Width(T1.xsize());
	BROCCOLI.SetT1Height(T1.ysize());
	BROCCOLI.SetT1Depth(T1.zsize());
	BROCCOLI.SetT1Timepoints(T1.tsize());
	BROCCOLI.SetT1VoxelSizeX(T1.xdim());
	BROCCOLI.SetT1VoxelSizeY(T1.ydim());
	BROCCOLI.SetT1VoxelSizeZ(T1.zdim());

	BROCCOLI.SetMNIWidth(mniVolume.xsize());
	BROCCOLI.SetMNIHeight(mniVolume.ysize());
	BROCCOLI.SetMNIDepth(mniVolume.zsize());
	BROCCOLI.SetMNIVoxelSizeX(mniVolume.xdim());
	BROCCOLI.SetMNIVoxelSizeY(mniVolume.ydim());
	BROCCOLI.SetMNIVoxelSizeZ(mniVolume.zdim());

	BROCCOLI.SetInterpolationMode(LINEAR);
	BROCCOLI.SetNumberOfIterationsForLinearImageRegistration(NUMBER_OF_ITERATIONS_FOR_LINEAR_IMAGE_REGISTRATION);
	BROCCOLI.SetNumberOfIterationsForNonLinearImageRegistration(0);
	BROCCOLI.SetImageRegistrationFilterSize(IMAGE_REGISTRATION_FILTER_SIZE);

	allocateLinearRegistrationFilter(filter);
	prepareLinearRegistrationFilters(filter);
	BROCCOLI.SetLinearImageRegistrationFilters(filter.h_Quadrature_Filter_1_Linear_Registration_Real, filter.h_Quadrature_Filter_1_Linear_Registration_Imag, filter.h_Quadrature_Filter_2_Linear_Registration_Real, filter.h_Quadrature_Filter_2_Linear_Registration_Imag, filter.h_Quadrature_Filter_3_Linear_Registration_Real, filter.h_Quadrature_Filter_3_Linear_Registration_Imag);

	BROCCOLI.SetCoarsestScaleT1MNI(COARSEST_SCALE);
	BROCCOLI.SetMMT1ZCUT(MM_T1_Z_CUT);

	BROCCOLI.SetTsigma(SIGMA);
	BROCCOLI.SetEsigma(SIGMA);
	BROCCOLI.SetDsigma(SIGMA);

	h_Interpolated_T1_Volume = (float *)malloc(mniVolume.xsize() * mniVolume.ysize() * mniVolume.zsize() * sizeof(float));
	h_Aligned_T1_Volume = (float *)malloc(mniVolume.xsize() * mniVolume.ysize() * mniVolume.zsize() * sizeof(float));

	BROCCOLI.SetOutputInterpolatedT1Volume(h_Interpolated_T1_Volume);
	BROCCOLI.SetOutputAlignedT1VolumeLinear(h_Aligned_T1_Volume);
	BROCCOLI.SetOutputT1MNIRegistrationParameters(h_Registration_Parameters);

	BROCCOLI.SetDoSkullstrip(MASK);
	BROCCOLI.SetDoSkullstripOriginal(MASK_ORIGINAL);

	BROCCOLI.SetSaveDisplacementField(WRITE_DISPLACEMENT_FIELD);
	BROCCOLI.SetSaveInterpolatedT1(WRITE_INTERPOLATED);
	BROCCOLI.SetSaveAlignedT1MNILinear(true);

	BROCCOLI.PerformRegistrationTwoVolumesWrapper();

	//printBufferErrors(BROCCOLI);
	//printRegistrationMatrix(h_Registration_Parameters);

	vector2volume(output, h_Aligned_T1_Volume);

	free(h_T1_Volume);
	free(h_MNI_Volume);

	free(h_Interpolated_T1_Volume);
	free(h_Aligned_T1_Volume);

	freeLinearRegistrationFilter(filter);
}

void BROCCOLIMotionCorrection(char *refFile, char *prefix, char *motionparms, int numFMRIs)
{
	BROCCOLI_LIB BROCCOLI(OPENCL_PLATFORM, OPENCL_DEVICE, 2, VERBOS); // 2 = Bash wrapper
	if (BROCCOLI.GetOpenCLInitiated())
	{
		float *h_Reference_Volume, *h_Motion_Parameters, *fmriVolume;
		broccoliLinearRegistrationFilters filter;
		volume<float> ref;

		load_volume(ref, string(refFile));
		h_Motion_Parameters = (float *)malloc(12 * sizeof(float));
		prepareBROCCOLIMotionCorrection(BROCCOLI, ref, h_Reference_Volume, h_Motion_Parameters, fmriVolume, filter);


		std::ofstream motion;
		motion.open(motionparms);

		if (motion.good())
			motion.precision(6);

		double startTime = GetWallTime();

		for (int t = 1; t <= numFMRIs; t++)
		{
			char filename[500], outputName[500];
			volume<float> epi;

			sprintf(filename, "%s%.5d.nii", prefix, t);
			//fprintf(stderr, "Trying to read %s \n", filename);
			load_volume(epi, string(filename));

			copyVolume2vector(epi, fmriVolume);
			BROCCOLI.PerformMotionCorrectionWrapper();
			vector2volume(epi, fmriVolume);

			sprintf(outputName, "%s%.5d_mc.nii", prefix, t);
			save_volume(epi, string(outputName));

			if (motion.good())
				motion << h_Motion_Parameters[4] << std::setw(2) << " " << -h_Motion_Parameters[3] << std::setw(2) << " " << h_Motion_Parameters[5] << std::setw(2) << " " << -h_Motion_Parameters[2] << std::setw(2) << " " << -h_Motion_Parameters[0] << std::setw(2) << " " << -h_Motion_Parameters[1] << std::endl;
		}
		double endTime = GetWallTime();
		printf("\nIt took %f seconds to perform just the motion correction\n", (float)(endTime - startTime));

		free(h_Motion_Parameters);
		free(h_Reference_Volume);
		free(fmriVolume);
		freeLinearRegistrationFilter(filter);
		motion.close();
	}
}


void prepareBROCCOLIMotionCorrection(BROCCOLI_LIB &BROCCOLI, volume<float> &ref, float *&h_Reference_Volume, float *&h_Motion_Parameters, float *&fmriVolume, broccoliLinearRegistrationFilters &filter)
{
	int64_t allocatedHostMemory = 0;

	allocatedHostMemory += volume2vector(ref, h_Reference_Volume);

	fmriVolume = (float *) malloc(allocatedHostMemory);
	// Set all necessary pointers and values
	BROCCOLI.SetInputfMRIVolumes(fmriVolume);

	BROCCOLI.SetAllocatedHostMemory(allocatedHostMemory);

	BROCCOLI.SetEPIWidth(ref.xsize());
	BROCCOLI.SetEPIHeight(ref.ysize());
	BROCCOLI.SetEPIDepth(ref.zsize());
	BROCCOLI.SetEPITimepoints(1);

	BROCCOLI.SetEPIVoxelSizeX(ref.xdim());
	BROCCOLI.SetEPIVoxelSizeY(ref.ydim());
	BROCCOLI.SetEPIVoxelSizeZ(ref.zdim());

	BROCCOLI.SetChangeMotionCorrectionReferenceVolume(true);
	BROCCOLI.SetMotionCorrectionReferenceVolume(h_Reference_Volume);

	BROCCOLI.SetImageRegistrationFilterSize(MOTION_CORRECTION_FILTER_SIZE);

	allocateLinearRegistrationFilter(filter);
	prepareLinearRegistrationFilters(filter);

	BROCCOLI.SetLinearImageRegistrationFilters(filter.h_Quadrature_Filter_1_Linear_Registration_Real, filter.h_Quadrature_Filter_1_Linear_Registration_Imag, filter.h_Quadrature_Filter_2_Linear_Registration_Real, filter.h_Quadrature_Filter_2_Linear_Registration_Imag, filter.h_Quadrature_Filter_3_Linear_Registration_Real, filter.h_Quadrature_Filter_3_Linear_Registration_Imag);
	BROCCOLI.SetNumberOfIterationsForMotionCorrection(NUMBER_OF_ITERATIONS_FOR_MOTION_CORRECTION);

	BROCCOLI.SetOutputMotionParameters(h_Motion_Parameters);
}

void printRegistrationMatrix(float *h_Registration_Parameters)
{
	printf("\nAffine registration parameters\n");
	printf(" %f %f %f %f\n", h_Registration_Parameters[3] + 1.0f, h_Registration_Parameters[4], h_Registration_Parameters[5], h_Registration_Parameters[0]);
	printf(" %f %f %f %f\n", h_Registration_Parameters[6], h_Registration_Parameters[7] + 1.0f, h_Registration_Parameters[8], h_Registration_Parameters[1]);
	printf(" %f %f %f %f\n", h_Registration_Parameters[9], h_Registration_Parameters[10], h_Registration_Parameters[11] + 1.0f, h_Registration_Parameters[2]);
	printf(" %f %f %f %f\n\n", 0.0f, 0.0f, 0.0f, 1.0f);
}

void writeRegistrationMatrix(char *output, float *h_Registration_Parameters)
{
	// Add the provided filename extension to the original filename, before the dot

	const char* extension = "_affinematrix.txt";

	std::ofstream matrix;
	matrix.open(output);
	if (matrix.good())
	{
		matrix.precision(6);

		matrix << h_Registration_Parameters[3] + 1.0f << std::setw(2) << " " << h_Registration_Parameters[4] << std::setw(2) << " " << h_Registration_Parameters[5] << std::setw(2) << " " << h_Registration_Parameters[0] << std::endl;
		matrix << h_Registration_Parameters[6] << std::setw(2) << " " << h_Registration_Parameters[7] + 1.0f << std::setw(2) << " " << h_Registration_Parameters[8] << std::setw(2) << " " << h_Registration_Parameters[1] << std::endl;
		matrix << h_Registration_Parameters[9] << std::setw(2) << " " << h_Registration_Parameters[10] << std::setw(2) << " " << h_Registration_Parameters[11] + 1.0f << std::setw(2) << " " << h_Registration_Parameters[2] << std::endl;
		matrix << 0.0f << std::setw(2) << " " << 0.0f << std::setw(2) << " " << 0.0f << std::setw(2) << " " << 1.0f << std::endl;

		matrix.close();
	}
	else
	{
		printf("Could not open %s for writing!\n", output);
	}
}

void printBufferErrors(BROCCOLI_LIB &BROCCOLI)
{
	// Print create buffer errors
	int* createBufferErrors = BROCCOLI.GetOpenCLCreateBufferErrors();
	for (int i = 0; i < BROCCOLI.GetNumberOfOpenCLKernels(); i++)
	{
		if (createBufferErrors[i] != 0)
		{
			printf("Create buffer error %i is %s \n", i, BROCCOLI.GetOpenCLErrorMessage(createBufferErrors[i]));
		}
	}

	// Print create kernel errors
	
	int* createKernelErrors = BROCCOLI.GetOpenCLCreateKernelErrors();
	for (int i = 0; i < BROCCOLI.GetNumberOfOpenCLKernels(); i++)
	{
		if (createKernelErrors[i] != 0)
		{
			printf("Create kernel error for kernel '%s' is '%s' \n", BROCCOLI.GetOpenCLKernelName(i), BROCCOLI.GetOpenCLErrorMessage(createKernelErrors[i]));
		}
	}

	// Print run kernel errors
	int* runKernelErrors = BROCCOLI.GetOpenCLRunKernelErrors();
	for (int i = 0; i < BROCCOLI.GetNumberOfOpenCLKernels(); i++)
	{
		if (runKernelErrors[i] != 0)
		{
			printf("Run kernel error for kernel '%s' is '%s' \n", BROCCOLI.GetOpenCLKernelName(i), BROCCOLI.GetOpenCLErrorMessage(runKernelErrors[i]));
		}
	}
}

void allocateLinearRegistrationFilter(broccoliLinearRegistrationFilters &filter)
{
	filter.h_Quadrature_Filter_1_Linear_Registration_Real = (float *)malloc(IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE * sizeof(float));
	filter.h_Quadrature_Filter_1_Linear_Registration_Imag = (float *)malloc(IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE * sizeof(float));

	filter.h_Quadrature_Filter_2_Linear_Registration_Real = (float *)malloc(IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE * sizeof(float));
	filter.h_Quadrature_Filter_2_Linear_Registration_Imag = (float *)malloc(IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE * sizeof(float));

	filter.h_Quadrature_Filter_3_Linear_Registration_Real = (float *)malloc(IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE * sizeof(float));
	filter.h_Quadrature_Filter_3_Linear_Registration_Imag = (float *)malloc(IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE * sizeof(float));
}

void freeLinearRegistrationFilter(broccoliLinearRegistrationFilters &filter)
{
	free(filter.h_Quadrature_Filter_1_Linear_Registration_Real);
	free(filter.h_Quadrature_Filter_1_Linear_Registration_Imag);

	free(filter.h_Quadrature_Filter_2_Linear_Registration_Real);
	free(filter.h_Quadrature_Filter_2_Linear_Registration_Imag);

	free(filter.h_Quadrature_Filter_3_Linear_Registration_Real);
	free(filter.h_Quadrature_Filter_3_Linear_Registration_Imag);
}

void prepareLinearRegistrationFilters(broccoliLinearRegistrationFilters &filter)
{
	std::string filter1RealLinearPathAndName;
	std::string filter1ImagLinearPathAndName;
	std::string filter2RealLinearPathAndName;
	std::string filter2ImagLinearPathAndName;
	std::string filter3RealLinearPathAndName;
	std::string filter3ImagLinearPathAndName;

	filter1RealLinearPathAndName.append(getenv("BROCCOLI_DIR"));
	filter1ImagLinearPathAndName.append(getenv("BROCCOLI_DIR"));
	filter2RealLinearPathAndName.append(getenv("BROCCOLI_DIR"));
	filter2ImagLinearPathAndName.append(getenv("BROCCOLI_DIR"));
	filter3RealLinearPathAndName.append(getenv("BROCCOLI_DIR"));
	filter3ImagLinearPathAndName.append(getenv("BROCCOLI_DIR"));

	filter1RealLinearPathAndName.append("filters/filter1_real_linear_registration.bin");
	filter1ImagLinearPathAndName.append("filters/filter1_imag_linear_registration.bin");
	filter2RealLinearPathAndName.append("filters/filter2_real_linear_registration.bin");
	filter2ImagLinearPathAndName.append("filters/filter2_imag_linear_registration.bin");
	filter3RealLinearPathAndName.append("filters/filter3_real_linear_registration.bin");
	filter3ImagLinearPathAndName.append("filters/filter3_imag_linear_registration.bin");

	// Read quadrature filters for linear registration, three real valued and three imaginary valued
	ReadBinaryFile(filter.h_Quadrature_Filter_1_Linear_Registration_Real, IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE, filter1RealLinearPathAndName.c_str());
	ReadBinaryFile(filter.h_Quadrature_Filter_1_Linear_Registration_Imag, IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE, filter1ImagLinearPathAndName.c_str());

	ReadBinaryFile(filter.h_Quadrature_Filter_2_Linear_Registration_Real, IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE, filter2RealLinearPathAndName.c_str());
	ReadBinaryFile(filter.h_Quadrature_Filter_2_Linear_Registration_Imag, IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE, filter2ImagLinearPathAndName.c_str());

	ReadBinaryFile(filter.h_Quadrature_Filter_3_Linear_Registration_Real, IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE, filter3RealLinearPathAndName.c_str());
	ReadBinaryFile(filter.h_Quadrature_Filter_3_Linear_Registration_Imag, IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE*IMAGE_REGISTRATION_FILTER_SIZE, filter3ImagLinearPathAndName.c_str());
}

void BROCCOLIEngineInterface::setRFI(char *rfiFileName)
{
	strcpy(rfiName, rfiFileName);
	load_volume(rfiVolume, string(rfiFileName));

	allocatedHostMemory += volume2vector(rfiVolume, h_Reference_Volume);
	allocatedHostMemory += volume2vector(rfiVolume, epiVolume);
	allocatedHostMemory += (1 + 2) * volume2vector(rfiVolume, h_Certainty);
	for (int t = 0; t < rfiVolume.xsize() * rfiVolume.ysize() * rfiVolume.zsize(); t++) h_Certainty[t] = 1.0f;

	BROCCOLI->SetInputfMRIVolumes(epiVolume);

	BROCCOLI->SetEPIWidth(rfiVolume.xsize());
	BROCCOLI->SetEPIHeight(rfiVolume.ysize());
	BROCCOLI->SetEPIDepth(rfiVolume.zsize());
	BROCCOLI->SetEPITimepoints(1);

	BROCCOLI->SetEPIVoxelSizeX(rfiVolume.xdim());
	BROCCOLI->SetEPIVoxelSizeY(rfiVolume.ydim());
	BROCCOLI->SetEPIVoxelSizeZ(rfiVolume.zdim());
}

void BROCCOLIEngineInterface::preparePipeline(float fwhm)
{
	prepareMotionCorrection();
	prepareSmoothing(fwhm);
}

void BROCCOLIEngineInterface::prepareSmoothing(float fwhm)
{
	fwhmValue = fwhm;
	BROCCOLI->SetAutoMask(false);
	BROCCOLI->SetInputCertainty(h_Certainty);

	BROCCOLI->SetEPISmoothingAmount(fwhm);
	BROCCOLI->SetAllocatedHostMemory(allocatedHostMemory);
}

void BROCCOLIEngineInterface::prepareRegistrationFilters()
{
	BROCCOLI->SetImageRegistrationFilterSize(MOTION_CORRECTION_FILTER_SIZE);
	allocateLinearRegistrationFilter(filter);
	prepareLinearRegistrationFilters(filter);
	BROCCOLI->SetLinearImageRegistrationFilters(filter.h_Quadrature_Filter_1_Linear_Registration_Real, filter.h_Quadrature_Filter_1_Linear_Registration_Imag, filter.h_Quadrature_Filter_2_Linear_Registration_Real, filter.h_Quadrature_Filter_2_Linear_Registration_Imag, filter.h_Quadrature_Filter_3_Linear_Registration_Real, filter.h_Quadrature_Filter_3_Linear_Registration_Imag);
}

void BROCCOLIEngineInterface::prepareMotionCorrection()
{
	// Set all necessary pointers and values
	BROCCOLI->SetChangeMotionCorrectionReferenceVolume(true);
	BROCCOLI->SetMotionCorrectionReferenceVolume(h_Reference_Volume);
	prepareRegistrationFilters();
	BROCCOLI->SetNumberOfIterationsForMotionCorrection(NUMBER_OF_ITERATIONS_FOR_MOTION_CORRECTION);

	BROCCOLI->SetOutputMotionParameters(h_Motion_Parameters);
}

void BROCCOLIEngineInterface::runSmoothing(volume<float> &epi)
{
	copyVolume2vector(epi, epiVolume);
	BROCCOLI->PerformSmoothingNormalizedHostWrapper();
	vector2volume(epi, epiVolume);
}

void BROCCOLIEngineInterface::createBROCCOLIObject(int platform, int device)
{
	BROCCOLI = new BROCCOLI_LIB(platform, device, 2, VERBOS);
}

void BROCCOLIEngineInterface::deallocateBROCCOLIObject()
{
	delete BROCCOLI;
	BROCCOLI = NULL;
}

void BROCCOLIEngineInterface::runMotionCorrection(volume<float> &epi)
{
	copyVolume2vector(epi, epiVolume);
	BROCCOLI->PerformMotionCorrectionWrapper();
	vector2volume(epi, epiVolume);
}

void BROCCOLIEngineInterface::pipelineEngine(char *epiFileName, char *outputFileName)
{
	volume <float>epi;
	load_volume(epi, string(epiFileName));
	pipelineEngine(epi);
	save_volume(epi, string(outputFileName));
}

void BROCCOLIEngineInterface::pipelineNormalEngine(char *epiFileName, char *outputFileName)
{
	char CmdLn[500], mcflirtParams[500];

	strcpy(mcflirtParams, "-plots -mats -rmsabs");

	sprintf(CmdLn, "mcflirt -in %s -reffile %s -out %s %s", epiFileName, rfiName, outputFileName, mcflirtParams);
	mcflirt(CmdLn);

	sprintf(CmdLn, "fslmaths %s -kernel gauss %f -fmean %s", outputFileName, fwhmValue / 2.3548, outputFileName);
	fslmaths(CmdLn);
}

void BROCCOLIEngineInterface::pipelineEngine(volume<float> &epi)
{
	// motion Correction
	runMotionCorrection(epi);
	// Baseline Calculation
	// Subtraction Baseline
	// Gaussian Smoothing
	runSmoothing(epi);
}