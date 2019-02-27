#include "fMRIUtils.h"

void znormalization(Matrix& mat, const int dim)
{
	Matrix matmean = mean(mat, dim);
	Matrix matstdev = stdev(mat, dim);

	if (dim == 1){
		for (int mr = 1; mr <= mat.Nrows(); mr++)
			mat.Row(mr) -= matmean.AsRow();

		for (int mc = 1; mc <= mat.Ncols(); mc++)
			mat.Column(mc) /= matstdev(1, mc);
	}
	else{
		for (int mc = 1; mc <= mat.Ncols(); mc++)
			mat.Column(mc) -= matmean.AsColumn();

		for (int mr = 1; mr <= mat.Nrows(); mr++)
			mat.Row(mr) /= matstdev(mr, 1);
	}
}

short getDataType(char *iname)
{
	string basename = fslbasename(string(iname));

	FSLIO* IP1;
	IP1 = FslOpen(basename.c_str(), "rb");
	if (IP1 == NULL) {
		cerr << "Cannot open volume " << basename << " for reading!\n";
		return -1;
	}

	short dtype;
	FslGetDataType(IP1, &dtype);
	float slope, intercept;
	int doscaling;
	doscaling = FslGetIntensityScaling(IP1, &slope, &intercept);
	if (doscaling == 1) { dtype = DT_FLOAT; }

	FslClose(IP1);
	FslFree(IP1);

	return NEWIMAGE::closestTemplatedType(dtype);
}

void volumeProperties(char *iname, int &runSize, short &datatype)
{
	string Onefile = iname;
	string basename = fslbasename(Onefile);

	FSLIO* IP1;
	IP1 = FslOpen(basename.c_str(), "rb");
	if (IP1 == NULL) {
		cerr << "Cannot open volume " << basename << " for reading!\n";
		return;
	}

	short dtype;
	FslGetDataType(IP1, &dtype);
	float slope, intercept;
	int doscaling;
	doscaling = FslGetIntensityScaling(IP1, &slope, &intercept);
	if (doscaling == 1) { dtype = DT_FLOAT; }
	runSize = IP1->niftiptr->nt;

	FslClose(IP1);
	FslFree(IP1);

	datatype = NEWIMAGE::closestTemplatedType(dtype);
}

template <class T>
int middlevol(char *iname, int middleIdx, char *oname)
{
	volume4D<T>vol;
	volume<T>outVol;

	read_volume4D(vol, string(iname));
	outVol.copyproperties(vol[middleIdx]);
	outVol.reinitialize(vol.xsize(), vol.ysize(), vol.zsize());
	outVol = vol[middleIdx];
	save_volume(outVol, oname);

	return 0;
}

int extractVolume(char *iname, char *oname, int idx, short datatype)
{
	if (datatype == DT_UNSIGNED_CHAR) return middlevol<char>(iname, idx, oname);
	else if (datatype == DT_SIGNED_SHORT) return middlevol<short>(iname, idx, oname);
	else if (datatype == DT_SIGNED_INT) return middlevol<int>(iname, idx, oname);
	else if (datatype == DT_FLOAT)  return middlevol<float>(iname, idx, oname);
	else if (datatype == DT_DOUBLE) return middlevol<double>(iname, idx, oname);
	return -1;
}

void splitVolume(char *fname, char *oname)
{
	volume4D<float>vol4D;
	volume<float>vol;

	read_volume4D(vol4D, string(fname));
	char saida[2048];


	FslSetOverrideOutputType(FSL_TYPE_NIFTI);
	for (int t = 0; t < vol4D.tsize(); t++)
	{
		vol = vol4D[t];
		sprintf(saida, "%s%.5d", oname, t + 1);
		save_volume(vol, string(saida));
	}
}
