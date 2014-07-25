#include "libvis/miscplot.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/t2z.h"
#include "libprob.h"
#include "parser.h"
#include "fslio/fslio.h"
#include "intervalos.h"

#include <cstring>

using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

extern "C" __declspec(dllexport) int _stdcall fslmaths(char *CmdLn);
extern "C" __declspec(dllexport) int _stdcall fslroi(char *CmdLn);
extern "C" __declspec(dllexport) int _stdcall susan(char *CmdLn);
template <class T> int _axial(char *innam, char *outnam);

int nifti_orientation_int( int ii )
{
   switch( ii ){
     case NIFTI_L2R: return 1 ;
     case NIFTI_R2L: return -1 ;
     case NIFTI_P2A: return 2 ;
     case NIFTI_A2P: return -2 ;
     case NIFTI_I2S: return 3 ;
     case NIFTI_S2I: return -3 ;
   }
   return 0;
}

void retornaijk(FSLIO* fslio, int *is, int *js, int *ks)
{
  int i, j, k;
  if (fslio!=NULL)
  {
     if( fslio->niftiptr->qform_code > 0 ) nifti_mat44_to_orientation( fslio->niftiptr->qto_xyz , &i,&j,&k );
     else if( fslio->niftiptr->sform_code > 0 ) nifti_mat44_to_orientation( fslio->niftiptr->sto_xyz , &i,&j,&k );
	 else 
	 {
		 i=NIFTI_R2L;
		 j=NIFTI_P2A;
		 k=NIFTI_I2S;

	 }
	 *is = nifti_orientation_int(i);
     *js = nifti_orientation_int(j);
     *ks = nifti_orientation_int(k);
  }
}

void retornaijk(char *innam, int *is, int *ys, int *zs)
{
  FSLIO* fslio=NULL;

  fslio = FslOpen(FslMakeBaseName(innam),"rb");
  retornaijk(fslio, is, ys, zs);
  if (fslio!=NULL) 
  {
	  FslClose(fslio);
	  FslFree(fslio);
  }

}

extern "C" __declspec(dllexport) void _stdcall CalculaBaseline4D(volume4D<float> *vol4D, int ini, int fim, volume<float> *&saida)
{
	saida = new volume<float>;
	saida->reinitialize(vol4D->xsize(), vol4D->ysize(), vol4D->zsize());
	saida->copyproperties((*vol4D)[0]);
	int total = fim-ini+1;
	
    for(int z=saida->minz(); z<=saida->maxz();z++){
       for(int y=saida->miny(); y<=saida->maxy();y++){
          for(int x=saida->minx(); x<=saida->maxx();x++)
	      {
			  float media=0;
			  for (int t=ini; t<=fim; t++)
			  {
				  media += (*vol4D)[t-1](x,y,z);
			  }
			  (*saida)(x,y,z) = media / total;
		  }
	   }
	}
}

template <class T>
int _axial2(char *innam, char *outnam)
{
  int newx, newy, newz;
  int i, j, k;
  int retval=0;
  char *basenam = FslMakeBaseName(innam);

  string inname=innam;
  string outname=outnam;

  volume4D<T> invol;
  read_orig_volume4D(invol,inname);
  FSLIO* fslio=NULL;

  fslio = FslOpen(basenam,"rb");
  free(basenam);
  if (fslio!=NULL)
  {
     Matrix affmat;
     retornaijk(fslio, &i, &j, &k);
	 FslClose(fslio);
	 FslFree(fslio);

	 if (i==1) newx=1;
	 else if (i==-1) newx=-1;
	 else if (j==-1) newx=-2;
	 else if (j==1) newx=2;
	 else if (k==-1) newx=-3;
	 else if (k==1) newx=3;

	 if (i==2) newy=1;
	 else if (i==-2) newy=-1;
	 else if (j==-2) newy=-2;
	 else if (j==2) newy=2;
	 else if (k==-2) newy=-3;
	 else if (k==2) newy=3;

	 if (i==3) newz=1;
	 else if (i==-3) newz=-1;
	 else if (j==-3) newz=-2;
	 else if (j==3) newz=2;
	 else if (k==-3) newz=-3;
	 else if (k==3) newz=3;

	 affmat = invol.swapmat(newx,newy,newz);
     for (int t=0; t<invol.tsize(); t++) {
       invol[t].basic_swapdimensions(newx,newy,newz,false);
     }

     if (outname!="") {
        retval = save_orig_volume4D(invol,outname);
     }
  }
  return retval;
}
/*
template <class T>
int _axial3(char *innam, char *outnam)
{
  int newx, newy, newz;
  int i, j, k;
  int retval=0;

  string inname=innam;
  string outname=outnam;

  volume4D<T> invol;
  read_orig_volume4D(invol,inname);
  FSLIO* fslio=NULL;

  fslio = FslOpen(FslMakeBaseName(innam),"rb");
  if (fslio!=NULL)
  {
     Matrix affmat;
     retornaijk(fslio, &i, &j, &k);
	 FslClose(fslio);
	 FslFree(fslio);

	 if (i==1) newx=1;
	 else if (i==-1) newx=-1;
	 else if (j==-1) newx=-2;
	 else if (j==1) newx=2;
	 else if (k==-1) newx=-3;
	 else if (k==1) newx=3;

	 if (i==2) newy=1;
	 else if (i==-2) newy=-1;
	 else if (j==-2) newy=-2;
	 else if (j==2) newy=2;
	 else if (k==-2) newy=-3;
	 else if (k==2) newy=3;

	 if (i==3) newz=1;
	 else if (i==-3) newz=-1;
	 else if (j==-3) newz=-2;
	 else if (j==3) newz=2;
	 else if (k==-3) newz=-3;
	 else if (k==3) newz=3;

	 affmat = invol.swapmat(newx,newy,newz);
     for (int t=0; t<invol.tsize(); t++) {
       invol[t].swapdimensions(newx,newy,newz);
     }

     if (outname!="") {
        retval = save_orig_volume4D(invol,outname);
     }
  }
  return retval;
}
*/

template <class T>
int _axial3(char *innam, char *outnam)
{
  int newx, newy, newz;
  int i, j, k;
  int retval=0;

  string inname=innam;
  string outname=outnam;

  volume4D<T> invol;
  read_orig_volume4D(invol,inname);
  FSLIO* fslio=NULL;

  fslio = FslOpen(FslMakeBaseName(innam),"rb");
  if (fslio!=NULL)
  {
     Matrix affmat;
     retornaijk(fslio, &i, &j, &k);
	 FslClose(fslio);
	 FslFree(fslio);

	 if (i==1) newx=1;
	 else if (i==-1) newx=-1;
	 else if (j==-1) newx=-2;
	 else if (j==1) newx=2;
	 else if (k==-1) newx=-3;
	 else if (k==1) newx=3;

	 if (i==2) newy=1;
	 else if (i==-2) newy=-1;
	 else if (j==-2) newy=-2;
	 else if (j==2) newy=2;
	 else if (k==-2) newy=-3;
	 else if (k==2) newy=3;

	 if (i==3) newz=1;
	 else if (i==-3) newz=-1;
	 else if (j==-3) newz=-2;
	 else if (j==3) newz=2;
	 else if (k==-3) newz=-3;
	 else if (k==3) newz=3;

	 affmat = invol.swapmat(newx,newy,newz);
     for (int t=0; t<invol.tsize(); t++) {
       invol[t].basic_swapdimensions(newx,newy,newz,false);
     }

     if (outname!="") {
        retval = save_orig_volume4D(invol,outname);
     }
  }
  return retval;
}

extern "C" __declspec(dllexport) void _stdcall brainthreshold(char *func, char* mascara, char *temp, char *saida, float brain_thres)
{
	float p2, p98;
	float thres;
	char cmd[500];
	char dquoted[10];
	volume4D<float> Vol;
	
	sprintf(dquoted, "\"");
	read_volume4D(Vol, string(func));
	p2 = Vol.percentile(0.02);
	p98 = Vol.percentile(0.98);
    thres = p2 + (brain_thres / 100) * (p98 - p2);
    sprintf(cmd, "fslmaths %s%s%s -thr %f -Tmin -bin %s%s%s -odt char", dquoted, func, dquoted, thres, dquoted, temp, dquoted);
	fslmaths(cmd);

//    sprintf(cmd, "fslmaths %s%s%s -dilF %s%s%s", dquoted, temp, dquoted, dquoted, temp, dquoted);
//	fslmaths(cmd);

    sprintf(cmd, "fslmaths %s%s%s -mas %s%s%s %s%s%s", dquoted, mascara, dquoted, dquoted, temp, dquoted, dquoted, saida, dquoted);
	fslmaths(cmd);
}

extern "C" __declspec(dllexport) float _stdcall percentile(char *vol, char* mascara, float perc)
{
	volume4D<float> Vol;
	volume<float> Mascara;

	read_volume4D(Vol, string(vol));
	if (mascara != NULL) 
	{
		read_volume(Mascara, string(mascara));
		return Vol.percentile(perc, Mascara);
	}
	else return Vol.percentile(perc);
}

extern "C" __declspec(dllexport) int _stdcall fazsusan(char *vol, char *saida, char *mascara, float FWHM)
{
	double p2, p98, p90, p50, threshold;
	volume<float> Vol;
	volume<float> Mascara;

	read_volume(Vol, string(vol));
	read_volume(Mascara, string(mascara));
	p2=Vol.percentile(0.02, Mascara);
	p50=Vol.percentile(0.50, Mascara);
	p90=Vol.percentile(0.90, Mascara);

	if ((p90-p2) != 0)
	{
	   char cmd[5000];
	   threshold = (p90-p2) * 0.75;
	   sprintf(cmd, "susan %s %f %f 3 1 1 %s %f %s", vol, threshold, FWHM/2.355, mascara, threshold, saida);
	   return susan(cmd);
	}
	else return -1;
}

extern "C" __declspec(dllexport) void _stdcall stdone(char *arq, char *mascara, char *saida)
{
	volume4D<float> tmpdata;
	read_volume4D(tmpdata,string(arq));
    volume<float> mask;   
	read_volume(mask,string(mascara));
    mask.binarise(0,mask.max()+1,exclusive); 

	Matrix data;
	data = tmpdata.matrix(mask);
	data =  SP(data,ones(data.Nrows(),1)*pow(stdev(data,1),-1));
	volume4D<float> tempVol;
	tempVol.setmatrix(data,mask);
	save_volume4D(tempVol,string(saida));
}

extern "C" __declspec(dllexport) void _stdcall demean(char *arq, char *saida)
{
   double mean;
   volume4D<float> input_volume;
   read_volume4D(input_volume,string(arq));
   volume4D<float> mask(input_volume);   
   mask.binarise(0,mask.max()+1,exclusive); 
   mean=input_volume.mean(mask);
   for(int t=0;t<input_volume.tsize();t++)     
      for(int z=0;z<input_volume.zsize();z++)
         for(int y=0;y<input_volume.ysize();y++)	    
	        for(int x=0;x<input_volume.xsize();x++)
               input_volume.value(x,y,z,t)=(float)(input_volume.value(x,y,z,t)-mean);
   save_volume4D(input_volume, string(saida));
}

template <class T>
int _samefov(char *base, char *mask, char *saida)
{
  short maskorig[5];
  short baseorig[5];
  int mx, my, mz;
  Matrix qmat, qbmat;
  volume<T> Base;
  read_volume(Base,string(base));

  volume<T> Mask;
  read_volume(Mask,string(mask));
 

   volume<T>masksamefov(Base.xsize(),Base.ysize(),Base.zsize());
   masksamefov.copyproperties(Base);

   // Volume Base
   if (Base.sform_code()!=NIFTI_XFORM_UNKNOWN) {
      qbmat = Base.sform_mat();
      baseorig[0]=abs(qbmat(1,4)) + 1;
      baseorig[1]=abs(qbmat(2,4)) + 1;
      baseorig[2]=abs(qbmat(3,4)) + 1;
   }
   else
   if (Base.qform_code()!=NIFTI_XFORM_UNKNOWN) {
      qbmat = Base.qform_mat();
      baseorig[0]=abs(qbmat(1,4)) + 1;
      baseorig[1]=abs(qbmat(2,4)) + 1;
      baseorig[2]=abs(qbmat(3,4)) + 1;
   };

   // Volume a Transformar
    if (Mask.sform_code()!=NIFTI_XFORM_UNKNOWN) {
      qmat = Mask.sform_mat();
	  maskorig[0]=abs(qmat(1,4)) + 1;
	  maskorig[1]=abs(qmat(2,4)) + 1;
	  maskorig[2]=abs(qmat(3,4)) + 1;
   }
   else
   if (Mask.qform_code()!=NIFTI_XFORM_UNKNOWN) {
      qmat = Mask.qform_mat();
	  maskorig[0]=abs(qmat(1,4)) + 1;
	  maskorig[1]=abs(qmat(2,4)) + 1;
	  maskorig[2]=abs(qmat(3,4)) + 1;
   }

   Matrix transf = qmat.i() * qbmat;
   float a11=transf(1,1), a12=transf(1,2), a13=transf(1,3), a14=transf(1,4),
	     a21=transf(2,1), a22=transf(2,2), a23=transf(2,3), a24=transf(2,4),
	     a31=transf(3,1), a32=transf(3,2), a33=transf(3,3), a34=transf(3,4);

   for(int z=Base.minz(); z<=Base.maxz();z++){
     for(int y=Base.miny(); y<=Base.maxy();y++){
       for(int x=Base.minx(); x<=Base.maxx();x++)
	   {
		   mx = x+(maskorig[0]-baseorig[0]);
		   my = y+(maskorig[1]-baseorig[1]);
		   mz = z+(maskorig[2]-baseorig[2]);

		   mx = a11*x + a12*y + a13*z + a14;
		   my = a21*x + a22*y + a23*z + a24;
		   mz = a31*x + a32*y + a33*z + a34;

		   if ((mx >= Mask.minx()) && (mx <= Mask.maxx()) && 
			   (my >= Mask.miny()) && (my <= Mask.maxy()) && 
			   (mz >= Mask.minz()) && (mz < Mask.maxz()))
			        masksamefov(x,y,z) = Mask(mx, my, mz);
		   else masksamefov(x,y,z) = 0;

	   }
	 }
   }
   save_volume(masksamefov, string(saida));
   return 1;
}

extern "C" __declspec(dllexport) int _stdcall samefov(char *base, char *mask, char *saida)
{
  string inname=mask;
  string basename = fslbasename(inname);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
	 FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return _samefov<char>(base, mask, saida);
     else if ( datatype==DT_SIGNED_SHORT ) return _samefov<short>(base, mask, saida);
     else if ( datatype==DT_SIGNED_INT ) return _samefov<int>(base, mask, saida);
     else if ( datatype==DT_FLOAT )  return _samefov<float>(base, mask, saida);
     else if ( datatype==DT_DOUBLE ) return _samefov<double>(base, mask, saida);
  }
  return -1;
}

template <class T>
int standarizavolume(char *base, char *mask, char *saida, int NN)
{
	char basetemp[500];
	char *basenam = FslMakeBaseName(base);
	volume<T> Base;
	read_volume(Base,string(base));
	sprintf(basetemp, "%st.nii", basenam);
	if (NN == 0)
	{
  	   resample<float>(mask, Base.xdim(), Base.ydim(), Base.zdim(), NN, saida);
	   _axial2<float>(saida, saida);
	   _axial2<float>(base, basetemp);
	   _samefov<float>(basetemp, saida, saida);
	}
	else
	{
  	   resample<T>(mask, Base.xdim(), Base.ydim(), Base.zdim(), NN, saida);
	   _axial2<T>(saida, saida);
	   _axial2<T>(base, basetemp);
	   _samefov<T>(basetemp, saida, saida);
	}
	remove(basetemp);
	free(basenam);
	return 1;
}

extern "C" __declspec(dllexport) int _stdcall StandarizaVolume(char *base, char *mask, char *saida, int NN)
{
  string inname=mask;
  string basename = fslbasename(inname);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return standarizavolume<char>(base, mask, saida, NN);
     else if ( datatype==DT_SIGNED_SHORT ) return standarizavolume<short>(base, mask, saida, NN);
     else if ( datatype==DT_SIGNED_INT ) return standarizavolume<int>(base, mask, saida, NN);
     else if ( datatype==DT_FLOAT )  return standarizavolume<float>(base, mask, saida, NN);
     else if ( datatype==DT_DOUBLE ) return standarizavolume<double>(base, mask, saida, NN);
  }
  return -1;
}

template <class T>
int aplicatransf3d(char *prefixovol, char *sufixovol, int inicio, int fim, char *volref, char *transmatrix, char *prefixo, int padding, int *x, int *y, int *z, int size, float FWHM)
{
   FILE *f;
   FWHM = (float) FWHM / (float) 2.3548;
   char arquivo[500];
   string inname;
   string refname=volref;
   char Saidatxt[500];
   sprintf(Saidatxt, "%s.txt", prefixo);

   f=fopen(Saidatxt, "wt+");
   if (f== NULL) return 100;

   volume<T> refvol, dummy;
   volume4D<T> invol;
   read_volume(refvol,refname);

   Matrix affmat(4,4);
   string matname = transmatrix;
   affmat = read_ascii_matrix(matname);

   volume<float> kernel;
   float xdim=0,ydim=0,zdim=0;
   dummy=refvol;
   for (int t=inicio; t<=fim; t++) 
   {
	   sprintf(arquivo, "%s%.5d%s", prefixovol, t, sufixovol);
	   inname=arquivo;
	   read_volume4D(invol, inname);

	   if (t==inicio)
	   {
		  kernel.destroy();
          kernel=gaussian_kernel3D(FWHM,invol[0].xdim(),invol[0].ydim(),invol[0].zdim());
	   }

	   invol[0].setpadvalue(invol[0].backgroundval());
	   invol.setextrapolationmethod(extraslice);
	   invol.setinterpolationmethod(trilinear);
	   invol=generic_convolve(invol,kernel,true,true);
       affine_transform(invol[0],dummy,affmat,padding);
	   dummy.setDisplayMaximumMinimum(0, 0);
	   for (int m=0; m<size;m++)
		   fprintf(f, "%f\n", (float) dummy(x[m], y[m], z[m]));
       if ((t % 20)==0) fflush(f);
   }
   fflush(f);
   fclose(f);
}

extern "C" __declspec(dllexport) int _stdcall AplicaTransformacao3D(char *prefixovol, char *sufixovol, int inicio, int fim, char *volref, char *transmatrix, char *prefixo, int padding, int *x, int *y, int *z, int size, float FWHM)
{
   char arquivo[500];
   sprintf(arquivo, "%s%.5d%s", prefixovol, inicio, sufixovol);
   string inname=arquivo;
   string basename = fslbasename(inname);

   FSLIO* IP1;
   IP1 = FslOpen(basename.c_str(),"rb");
   if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
   }
   else
   {
      short tipo;
      FslGetDataType(IP1,&tipo);
      FslClose(IP1);
      FslFree(IP1);

      short datatype=NEWIMAGE::closestTemplatedType(tipo);
      if ( datatype==DT_UNSIGNED_CHAR ) return aplicatransf3d<char>(prefixovol, sufixovol, inicio, fim, volref, transmatrix, prefixo, padding, x, y, z, size, FWHM);
      else if ( datatype==DT_SIGNED_SHORT ) return aplicatransf3d<short>(prefixovol, sufixovol, inicio, fim, volref, transmatrix, prefixo, padding, x, y, z, size, FWHM);
      else if ( datatype==DT_SIGNED_INT ) return aplicatransf3d<int>(prefixovol, sufixovol, inicio, fim, volref, transmatrix, prefixo, padding, x, y, z, size, FWHM);
      else if ( datatype==DT_FLOAT )  return aplicatransf3d<float>(prefixovol, sufixovol, inicio, fim, volref, transmatrix, prefixo, padding, x, y, z, size, FWHM);
      else if ( datatype==DT_DOUBLE ) return aplicatransf3d<double>(prefixovol, sufixovol, inicio, fim, volref, transmatrix, prefixo, padding, x, y, z, size, FWHM);
   }
   return -1;
}

template <class T>
int aplicatransf4d(char *vol4d, char *volref, char *transmatrix, char *prefixo, int padding, int *x, int *y, int *z, int size)
{
   FILE *f;
   string inname=vol4d;
   string refname=volref;
   char Saidatxt[500];
   sprintf(Saidatxt, "%s.txt", prefixo);

   f=fopen(Saidatxt, "wt+");
   if (f== NULL) return 100;

   volume4D<T> invol;
   volume<T> refvol, dummy;
   read_volume4D(invol,inname);
   for (int t=0; t<invol.tsize(); t++) {
      invol[t].setpadvalue(invol[t].backgroundval());
   }
   invol.setextrapolationmethod(extraslice);
   invol.setinterpolationmethod(trilinear);


    // old form used a volume number
    //    refvol = invol[atoi(refname.c_str())];  
    read_volume(refvol,refname);

    Matrix affmat(4,4);
    string matname = transmatrix;
    affmat = read_ascii_matrix(matname);
    dummy = refvol;
	for (int m=invol.mint(); m<=invol.maxt(); m++) 
	{
       //dummy = refvol;
       affine_transform(invol[m],dummy,affmat,padding);
	   dummy.setDisplayMaximumMinimum(0, 0);
	   char nomesaida[500];
	   sprintf(nomesaida, "%s%.5d.nii", prefixo, m+1);
	   string saida=nomesaida;
//	   fprintf(f, "\n");
//	   fprintf(f, "%d\n", m)+1;
	   for (int t=0; t<size;t++)
	   {
		   fprintf(f, "%f\n", dummy(x[t], y[t], z[t]));
	   };
//	   fprintf(f, "\n");
	   if ((m % 20)==0) fflush(f);
//	   save_volume(dummy, saida);
	}
	fflush(f);
	fclose(f);
}

extern "C" __declspec(dllexport) int _stdcall AplicaTransformacao4D(char *vol4d, char *volref, char *transmatrix, char *prefixo, int padding, int *x, int *y, int *z, int size)
{
  string inname=vol4d;
  string basename = fslbasename(inname);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return aplicatransf4d<char>(vol4d, volref, transmatrix, prefixo, padding, x, y, z, size);
     else if ( datatype==DT_SIGNED_SHORT ) return aplicatransf4d<short>(vol4d, volref, transmatrix, prefixo, padding, x, y, z, size);
     else if ( datatype==DT_SIGNED_INT ) return aplicatransf4d<int>(vol4d, volref, transmatrix, prefixo, padding, x, y, z, size);
     else if ( datatype==DT_FLOAT )  return aplicatransf4d<float>(vol4d, volref, transmatrix, prefixo, padding, x, y, z, size);
     else if ( datatype==DT_DOUBLE ) return aplicatransf4d<double>(vol4d, volref, transmatrix, prefixo, padding, x, y, z, size);
  }
  return -1;
}

ReturnMatrix do_convolve(const ColumnVector input, const ColumnVector kernel, const int phase, const int renorm)
{
  ColumnVector output(input);
  
  output=0;

  for(int t=0; t<input.Nrows(); t++)
    {
      float kernel_norm=0;
      for(int i=MAX(0,1+t+phase-input.Nrows()); i<MIN(kernel.Nrows(),t+phase+1); i++)
	{
	  output(t+1) += input(t+phase-i+1) * kernel(i+1);
	  kernel_norm += kernel(i+1);
	}
      if (renorm)
	output(t+1) /= kernel_norm;
    }

  output.Release();
  return output;
}

// }}}
// {{{ resample down in time 

/* sample in the MIDDLE of the upsampled period */

void do_resample(const ColumnVector input, Matrix &output, const int real_ev, float trmult, int negpts)
{
  for(int t=0;t<output.Nrows();t++)
    output(t+1,real_ev+1) = input(((int)((t+0.5)*trmult))+negpts+1);
}

ReturnMatrix mygammapdf(const int npts, const float mult, const float delay, const float sigma)
{
  ColumnVector grot(npts);
  
  for(int i=0; i<npts; i++)
    grot(i+1)=i/mult;

  grot = gammapdf(grot.t(),delay,sigma*sigma).t();

  grot.Release();
  return grot;
}

extern "C" __declspec(dllexport) float * _stdcall AlocaFloatArray(int tamanho)
{
	return (float *) malloc(tamanho * sizeof(float));
}

extern "C" __declspec(dllexport) void _stdcall LiberaFloatArray(float *tam)
{
	free(tam);
}

extern "C" __declspec(dllexport) int _stdcall convolui(int tamanho, float Mult, float *valores, float *saida)
{
#define DT                0.05 /* temporal sampling for initial model, in seconds */
#define NEGSECS           30.0 /* amount of negative modelling allowed, for custom 3 etc, in seconds */
    float mult = 1/DT;
	float sigma1=2.449, delay1=6, // first gamma
	sigma2=4, delay2=16,    // second gamma
	ratio=6;                // hrf = gammapdf1 - gammapdf2/ratio;
	float maxconvwin=0;
	float tr=2;
	float trmult=tr*mult;
	ColumnVector input;
	float convolve_phase=1;

	int fw = (int)((delay2 + sigma2*5));
	mult=Mult;
	maxconvwin=MAX(fw,tamanho);

	ColumnVector cX = mygammapdf(fw,mult,delay1,sigma1) - mygammapdf(fw,mult,delay2,sigma2)/ratio;
	input.ReSize(tamanho);
	for (int i=1;i<=input.Nrows();i++) input(i)=valores[i-1];

	ColumnVector oX=do_convolve(input,cX,convolve_phase,1);
	for (int i=1;i<=tamanho;i++) saida[i-1]=oX(i);
	
    //Matrix outp(tamanho, 1);
//	do_resample(oX, outp, 0, trmult, (int)(NEGSECS*mult));

//	for (int i=1;i<=tamanho;i++) saida[i-1]=outp(i, 1);

	return 0;
}

typedef struct
{
	int de, para;
} Mapeamento;

int LeMapeamentos(char *arquivo, Mapeamento **mapeamentos)
{
	ifstream in;
	string linha;
	in.open(arquivo);
	
    int tamanho=0;
	*mapeamentos=NULL;
	while (getline(in, linha))
	{
       *mapeamentos= (Mapeamento *) realloc(*mapeamentos, ++tamanho*sizeof(Intervalo));
	   if (*mapeamentos != NULL)
	   {
           char num[20];
		   int pos = linha.find("=");
		   int pos2= linha.length();
		   linha.copy(num, pos);
		   num[pos]='\0';
           (*mapeamentos)[tamanho-1].de = atoi(num);
		   linha.copy(num, pos2-pos-1, pos+1);
           (*mapeamentos)[tamanho-1].para = atoi(num);
	   }
	}
	in.close();
	return tamanho;
}

extern "C" __declspec(dllexport) int _stdcall Gera4DMedioParaVoxelEstavel(char *Arq4D, char *Saida, char *Design, char *ArqMap)
{
   volume4D<float> arq4d;
   if (Arq4D != NULL) 
   {
      string arq4D = Arq4D;
	  read_volume4D(arq4d, arq4D);
   }
   
   Mapeamento *mapeamentos;
   int tamMap=LeMapeamentos(ArqMap, &mapeamentos);

   Intervalo *intervalos;
   int tamInt=LeIntervalo(Design, &intervalos);
    

   volume4D<float>resultado(arq4d.xsize(),arq4d.ysize(),arq4d.zsize(),tamMap);
   resultado.copyproperties(arq4d);
   for(int i=0; i < tamMap; i++)
   {
	   int ini = intervalos[mapeamentos[i].para-1].inicio;
	   int fim = intervalos[mapeamentos[i].para-1].fim;
	   for (int j=ini; j<=fim; j++) resultado[i] += arq4d[j-1];
	   resultado[i] /= (float) (fim-ini+1);
   }
   string outname = Saida;
   save_volume4D(resultado, outname); 
   free(mapeamentos);
   free(intervalos);
   return 0;
}

extern "C" __declspec(dllexport) int _stdcall CalculaAtivacao(int Inicio, int ini, int fim, int tamb, char *prefixo, char *Basal, char *Mascara, char *Saida, int Width)
{
   char tmpname[5000];

   volume<float> basal;
   if (Basal != NULL) 
   {
      string bName = Basal;
      read_volume(basal, bName);
   }

   char formatstring[50];

   sprintf(formatstring, " s .%dd.nii", Width);
   formatstring[0] = '%';
   formatstring[2] = '%';

   volume<float> auxvol;
//   sprintf(tmpname, formatstring, prefixo, ini);
   sprintf(tmpname, prefixo, ini);
   string auxn = tmpname;
   read_volume(auxvol, auxn);

   volume4D<float>ativacao(auxvol.xsize(),auxvol.ysize(),auxvol.zsize(),fim-ini+1);
   ativacao.copyproperties(auxvol);
   for(int i=ini; i <= fim; i++)
   {
	   int comeco=Max(Inicio, i-tamb+1);
	   ativacao[i-ini] = 0;
	   for (int j=comeco; j<=i; j++)
	   {
  	      volume<float> tmp;
	      sprintf(tmpname, prefixo, j);
		  string Tmpname = tmpname;
          read_volume(tmp, Tmpname);
		  ativacao[i-ini] += tmp;
	   }
	   ativacao[i-ini] /= (float) (i-comeco+1);
	   if (Basal != NULL) ativacao[i-ini]-=basal;
   }
   if (Mascara != NULL)
   {
	   volume<float> mascara;
	   string mascarafile = Mascara;
	   read_volume(mascara, mascarafile);
	   mascara.binarise(0,mascara.max()+1,exclusive);
	   for (int i=ini; i<=fim; i++) ativacao[i-ini]*=mascara;
   }
   string outname = Saida;
   save_volume4D(ativacao, outname); 
   return 0;
}

extern "C" __declspec(dllexport) int _stdcall mem_AbreMascaraBinarizada(char *Mascara, volume<float> *&mascara)
{
	mascara = new volume<float>;
	string mask = Mascara;
	read_volume(*mascara, mask);
	mascara->binarise(0,mascara->max()+1,exclusive);
	return 0;
}

extern "C" __declspec(dllexport) int _stdcall mem_CalculaAtivacao4DVolume(int Inicio, int ini, int tamb, volume4D<float> *vol4D, volume<float> *basal, volume<float> *mascara, volume<float> *&ativacao)
{
	ativacao = new volume<float>;
	ativacao->reinitialize(vol4D->xsize(),vol4D->ysize(),vol4D->zsize());
    ativacao->copyproperties((*vol4D)[0]);
    int comeco=Max(Inicio, ini-tamb+1);
	*ativacao = 0;
	for (int j=comeco; j<=ini; j++) *ativacao += (*vol4D)[j-Inicio];
	*ativacao /= (float) (ini-comeco+1);
	if (basal != NULL) *ativacao-=*basal;

	if (mascara != NULL) (*ativacao)*=(*mascara);
	return 0;
}

extern "C" __declspec(dllexport) int _stdcall mem_CalculaAtivacao4D(int Inicio, int ini, int fim, int tamb, volume4D<float> *vol4D, char *Basal, char *Mascara, char *Saida)
{
   char tmpname[5000];
   int inc=Inicio;

   volume<float> basal;
   if (Basal != NULL) 
   {
      string bName = Basal;
      read_volume(basal, bName);
   }
   
   volume4D<float>ativacao(vol4D->xsize(),vol4D->ysize(),vol4D->zsize(),fim-ini+1);
   ativacao.copyproperties(*vol4D);
   for(int i=ini; i <= fim; i++)
   {
	   int comeco=Max(Inicio, i-tamb+1);
	   ativacao[i-ini] = 0;
	   for (int j=comeco; j<=i; j++) ativacao[i-ini] += (*vol4D)[j-inc];
	   ativacao[i-ini] /= (float) (i-comeco+1);
	   if (Basal != NULL) ativacao[i-ini]-=basal;
   }
   if (Mascara != NULL)
   {
	   volume<float> mascara;
	   string mascarafile = Mascara;
	   read_volume(mascara, mascarafile);
	   mascara.binarise(0,mascara.max()+1,exclusive);
	   for (int i=ini; i<=fim; i++) ativacao[i-ini]*=mascara;
   }
   string outname = Saida;
   save_volume4D(ativacao, outname); 
   return 0;
}

extern "C" __declspec(dllexport) int _stdcall CalculaAtivacao4D(int Inicio, int ini, int fim, int tamb, char *Vol4D, char *Basal, char *Mascara, char *Saida)
{
   volume<float> basal;
   if (Basal != NULL) 
   {
      string bName = Basal;
      read_volume(basal, bName);
   }
   
   volume4D<float> vol4D;
   if (Vol4D != NULL) 
   {
      string arqvol4D = Vol4D;
      read_volume4D(vol4D, arqvol4D);
   }
   return mem_CalculaAtivacao4D(Inicio, ini, fim, tamb, &vol4D, Basal, Mascara, Saida);
}

extern "C" __declspec(dllexport) int _stdcall Filtra4D(char *Vol4D, int *indices, int tamindices, char *Saida)
{
   volume4D<float> vol4D;
   if (Vol4D != NULL) 
   {
      string arqvol4D = Vol4D;
      read_volume4D(vol4D, arqvol4D);
   }

   volume4D<float>ativacao(vol4D.xsize(),vol4D.ysize(),vol4D.zsize(),tamindices);
   ativacao.copyproperties(vol4D);
   for(int i=0; i < tamindices; i++)  ativacao[i] = vol4D[indices[i]-1];
   string outname = Saida;
   save_volume4D(ativacao, outname); 
   return 0;
}

int pegaltimobaseline(int intervalo, Intervalo *intervalos, int tamintervalos, char *condicaobasal)
{
	int baseline = -1;
	if (condicaobasal != NULL)
	{
		if (intervalo>0)
		{
		   for (int t=intervalo-1;t>=0;t--)
		   {
			  if (strcmp(condicaobasal, intervalos[t].nome)==0)
			  {
				  baseline=t+1;
				  break;
			  }
		   }
		}
		if (baseline==-1)
		{
			for (int t=intervalo;t<tamintervalos;t++)
			{
			   if (strcmp(condicaobasal, intervalos[t].nome)==0)
			   {
				   baseline=t+1;
				   break;
			   }
			}
		}
	}
	else baseline = intervalo; //ultimo intervalo (os arquivos comecam de 1, diferente da variavel intervalo, que comeca de zero)
	return baseline;
}

extern "C" __declspec(dllexport) int _stdcall NormalizaPorBaseLine(char *Vol4D, char *arqintervalos, char *condicaobasal, char *Saida)
{
   Intervalo *intervalos;
   int tamintervalos=LeIntervalo(arqintervalos, &intervalos);
   volume<float> *basal=NULL, *newbasal=NULL;
   int numbasal=-1;
   volume4D<float> vol4D;
   if (Vol4D != NULL) 
   {
      string arqvol4D = Vol4D;
      read_volume4D(vol4D, arqvol4D);
   }
   int fim = vol4D.tsize();
   for(int i=0; i < fim; i++)
   {
	   int interv=retornaintervalo(i+1, intervalos, tamintervalos);
	   if (interv>=0)
	   {
		   if ((i+1)== intervalos[interv].inicio)
		   {
   	          if (strcmpi(intervalos[interv].nome, condicaobasal)==0)
		      {
			      if (newbasal !=NULL) delete newbasal;
			      CalculaBaseline4D(&vol4D, intervalos[interv].inicio, intervalos[interv].fim, newbasal);
			  }
		   }
		   if (basal != NULL) vol4D[i]-=*basal;
		   else vol4D[i]=0;

		   if ((i+1)== intervalos[interv].fim)
		   {
			   if (basal != NULL) delete basal;
			   basal=newbasal;
			   newbasal=NULL;
		   }

	   }
   }
   string outname = Saida;
   save_volume4D(vol4D, outname); 
   free(intervalos);
   if (newbasal !=NULL) delete newbasal;
   if (basal !=NULL) delete basal;

   return 0;
}

extern "C" __declspec(dllexport) int _stdcall retornaparams(char *matrixfile, char *volnam, float *rx, float *ry, float *rz, float *tx, float *ty, float *tz)
{
   Matrix affmat(4,4);
   ColumnVector params(12), cor(3);
   cor=0;  // centre of rotations
   affmat = read_ascii_matrix(matrixfile);
   if (affmat.Nrows()<4) return -2;
   if (volnam) 
   {
      string volname=volnam;
      volume<float> testvol;
      if (read_volume(testvol,volname)<0)  return -1;
      cor = testvol.cog("scaled_mm");
   };
   decompose_aff(params,affmat,cor,rotmat2euler);

   *rx=params(1);
   *ry=params(2);
   *rz=params(3);

   *tx=params(4);
   *ty=params(5);
   *tz=params(6);
   return 0;
}

template <class T>
int _axial(char *innam, char *outnam)
{
  int newx, newy, newz;
  int i, j, k;
  int retval=0;

  string inname=innam;
  string outname=outnam;

  volume4D<T> invol;
  read_orig_volume4D(invol,inname);
  FSLIO* fslio=NULL;

  fslio = FslOpen(FslMakeBaseName(innam),"rb");
  if (fslio!=NULL)
  {
     Matrix affmat;
     if( fslio->niftiptr->qform_code > 0 ) nifti_mat44_to_orientation( fslio->niftiptr->qto_xyz , &i,&j,&k );
     else if( fslio->niftiptr->sform_code > 0 ) nifti_mat44_to_orientation( fslio->niftiptr->sto_xyz , &i,&j,&k );
	 else 
	 {
		 i=NIFTI_R2L;
		 j=NIFTI_P2A;
		 k=NIFTI_I2S;

	 }
	 FslClose(fslio);
	 FslFree(fslio);

	 i = nifti_orientation_int(i);
     j = nifti_orientation_int(j);
     k = nifti_orientation_int(k);

	 if (i==1) newx=-1;
	 else if (i==-1) newx=1;
	 else if (j==-1) newx=2;
	 else if (j==1) newx=-2;
	 else if (k==-1) newx=3;
	 else if (k==1) newx=-3;

	 if (i==2) newy=1;
	 else if (i==-2) newy=-1;
	 else if (j==-2) newy=-2;
	 else if (j==2) newy=2;
	 else if (k==-2) newy=-3;
	 else if (k==2) newy=3;

	 if (i==3) newz=1;
	 else if (i==-3) newz=-1;
	 else if (j==-3) newz=-2;
	 else if (j==3) newz=2;
	 else if (k==-3) newz=-3;
	 else if (k==3) newz=3;

	 affmat = invol.swapmat(newx,newy,newz);

     if (affmat.Determinant()<0.0) {
        cout << "WARNING:: Flipping Left/Right orientation (as det < 0)" << endl;
     }
  
     for (int t=0; t<invol.tsize(); t++) {
       invol[t].basic_swapdimensions(newx,newy,newz,false);
     }

     if (outname!="") {
        retval = save_orig_volume4D(invol,outname);
     }
  }
  return retval;
}

extern "C" __declspec(dllexport) int _stdcall axial2(char *innam, char *outnam)
{
  string inname=innam;
  string basename = fslbasename(innam);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return _axial2<char>(innam, outnam);
     else if ( datatype==DT_SIGNED_SHORT ) return _axial2<short>(innam, outnam);
     else if ( datatype==DT_SIGNED_INT ) return _axial2<int>(innam, outnam);
     else if ( datatype==DT_FLOAT )  return _axial2<float>(innam, outnam);
     else if ( datatype==DT_DOUBLE ) return _axial2<double>(innam, outnam);
  }
  return -1;
}

extern "C" __declspec(dllexport) int _stdcall axial3(char *innam, char *outnam)
{
  string inname=innam;
  string basename = fslbasename(innam);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return _axial3<char>(innam, outnam);
     else if ( datatype==DT_SIGNED_SHORT ) return _axial3<short>(innam, outnam);
     else if ( datatype==DT_SIGNED_INT ) return _axial3<int>(innam, outnam);
     else if ( datatype==DT_FLOAT )  return _axial3<float>(innam, outnam);
     else if ( datatype==DT_DOUBLE ) return _axial3<double>(innam, outnam);
  }
  return -1;
}

extern "C" __declspec(dllexport) int _stdcall axial(char *innam, char *outnam)
{
  string inname=innam;
  string basename = fslbasename(innam);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return _axial<char>(innam, outnam);
     else if ( datatype==DT_SIGNED_SHORT ) return _axial<short>(innam, outnam);
     else if ( datatype==DT_SIGNED_INT ) return _axial<int>(innam, outnam);
     else if ( datatype==DT_FLOAT )  return _axial<float>(innam, outnam);
     else if ( datatype==DT_DOUBLE ) return _axial<double>(innam, outnam);
  }
  return -1;
}

extern "C" __declspec(dllexport) FSLIO * _stdcall fslioopen(char *arquivo)
{
  FSLIO *src = FslOpen(arquivo, "r");
  return src;  
}

extern "C" __declspec(dllexport) void _stdcall fslioclose(FSLIO *src)
{
  FslClose(src);
  FslFree(src);
}

template <class T>
int fmrib_main2(int argc,char *argv[], FSLIO *src)
{
  string newx=argv[2], newy=argv[3], newz=argv[4];
  string inname=argv[1];
  string outname="";

  if (argc>5) {
    outname=argv[5];
  };

  volume4D<T> invol;
  read_orig_volume4D(invol,inname);

  Matrix affmat;
  affmat = invol.swapmat(newx,newy,newz);


  invol.swapdimensions(newx,newy,newz);

  /*
  if (affmat.Determinant()<0.0) {
    invol.swapdimensions(-1,2,3);
  }
  */

  int retval=0;

  if (src!=NULL) {
      float vx, vy, vz, tr;
	  short scode, qcode;
      mat44 smat, qmat;
	  Matrix sMat(4, 4), qMat(4,4);

	  {
         scode = FslGetStdXform(src,&smat);
         qcode = FslGetRigidXform(src,&qmat);
         FslGetVoxDim(src, &vx, &vy, &vz, &tr);
		 for (int i=0; i<4; i++)
		 {
			 for (int j=0; j<4;j++)
			 {
		        sMat.element(i,j) = smat.m[i][j];
		        qMat.element(i,j) = qmat.m[i][j];
			 }
		 }

		 invol.set_qform(qcode, qMat);
		 invol.set_sform(scode, sMat);
		 invol.setxdim(vx);
		 invol.setydim(vy);
		 invol.setzdim(vz);
		 invol.setTR(tr);
	  }
  }

  if (outname!="") {
    retval = save_orig_volume4D(invol,outname);
  }
  return retval;
}

extern "C" __declspec(dllexport) void _stdcall mem_desalocamatriz(vector<Matrix> * mat)
{
	delete mat;
}

/*
extern "C" __declspec(dllexport) void _stdcall mem_desalocamatriz(Matrix * mat)
{
	delete mat;
}
*/

extern "C" __declspec(dllexport) void _stdcall mem_loadvolume(char *input, volume<float>* &vol)
{
	string inputname=input;
	vol = new volume<float>;
	load_volume(*vol, inputname);
}

extern "C" __declspec(dllexport) void _stdcall mem_loadvolume4D(char *input, volume4D<float>* &vol)
{
	string inputname=input;
	vol = NULL;
    // verificando se consigo abrir o arquivo
    FSLIO *IP1;
    IP1 = NewFslOpen(input, "r");
    if (IP1==0) return;
    else 
	{
		FslClose(IP1);
		FslFree(IP1);
	}

	vol = new volume4D<float>;
	load_volume4D(*vol, inputname);
}

extern "C" __declspec(dllexport) void _stdcall Zera4D(volume4D<float>* vol)
{
	for (int t=0; t<vol->tsize(); t++) (*vol)[t] = 0;
}

extern "C" __declspec(dllexport) void _stdcall mem_savevolume4D(char *Saida, volume4D<float>* volume)
{
	string outname=Saida;
	save_volume4D(*volume, outname);
}

extern "C" __declspec(dllexport) void _stdcall mem_desalocavolume4D(volume4D<float>* volume)
{
	delete volume;
}

extern "C" __declspec(dllexport) void _stdcall mem_desalocavolume(volume<float>* volume)
{
	delete volume;
}

extern "C" __declspec(dllexport) int _stdcall concatTransformations(char *volA, char *volB, char *affMatBC, char *affMatAB, char *outMatFile)
{
   volume<float> volumeA, volumeB;
   Matrix affmatAB(4,4), affmatBC(4,4);
   affmatAB = read_ascii_matrix(affMatAB);
   affmatBC = read_ascii_matrix(affMatBC);

   read_volume(volumeA, string(volA));
   read_volume(volumeB, string(volB));

   // iaffbig goes from output mm coords to input (reference) mm coords
   Matrix iaffbig = affmatAB;
   // check the left-right data orientations of the images and modify
   //   the transformation matrix to flip to radiological coords if necessary
   if (volumeA.left_right_order()==FSL_NEUROLOGICAL) 
   {
	  iaffbig = iaffbig * volumeA.swapmat(-1,2,3);
   }
   if (volumeB.left_right_order()==FSL_NEUROLOGICAL) 
   {
	  iaffbig = volumeB.swapmat(-1,2,3) * iaffbig;
   }
   
   // convert iaffbig to go from output voxel coords to input (reference) voxel coords
   iaffbig = volumeB.sampling_mat() * iaffbig * volumeA.sampling_mat().i();
   iaffbig = affmatBC * iaffbig;

   write_ascii_matrix(iaffbig, outMatFile);
   return 0;
}

extern "C" __declspec(dllexport) int _stdcall mem_fslswapdim_rt(char *CmdLn, FSLIO *src, volume4D<float>* &saida)
{
   int argc;
   char **argv;
  
   saida = NULL;
   parser(CmdLn, argc, argv);

   string newx=argv[2], newy=argv[3], newz=argv[4];
   string inname=argv[1];

   // verificando se consigo abrir o arquivo
   FSLIO *IP1;
   IP1 = NewFslOpen(argv[1], "r");
   if (IP1==0) {
	   freeparser(argc, argv);
	   return -1;
   }
   else 
   {
	   FslClose(IP1);
	   FslFree(IP1);
   }

   saida = new volume4D<float>;
   read_orig_volume4D(*saida,inname);

   Matrix affmat;
   affmat = saida->swapmat(newx,newy,newz);

   saida->swapdimensions(newx,newy,newz);
//   if (affmat.Determinant()<0.0)
//       saida->swapdimensions(-1,2,3);


   int retval=0;

   if (src!=NULL) 
   {
      float vx, vy, vz, tr;
	  short scode, qcode;
      mat44 smat, qmat;
	  Matrix sMat(4, 4), qMat(4,4);
      scode = FslGetStdXform(src,&smat);
      qcode = FslGetRigidXform(src,&qmat);
      FslGetVoxDim(src, &vx, &vy, &vz, &tr);
	  for (int i=0; i<4; i++)
	  {
		 for (int j=0; j<4;j++)
		 {
		    sMat.element(i,j) = smat.m[i][j];
		    qMat.element(i,j) = qmat.m[i][j];
		 }
	  }
	  saida->set_qform(qcode, qMat);
	  saida->set_sform(scode, sMat);
	  saida->setxdim(vx);
	  saida->setydim(vy);
	  saida->setzdim(vz);
	  saida->setTR(tr);
   }
   else 
   {
	   delete saida;
	   saida = NULL;
   }
   freeparser(argc, argv);
   return retval;
}

extern "C" __declspec(dllexport) int _stdcall fslswapdim_rt(char *CmdLn, FSLIO *src)
{
  int r;
  int argc;
  char **argv;
  
  parser(CmdLn, argc, argv);

  string progname=argv[0];
  if (argc<5) { 
    //print_usage(progname);
    r= -1; 
  }
  else
  {
    string inname = argv[1];
    short datatype1=NEWIMAGE::closestTemplatedType(dtype(inname));
    if ( datatype1==DT_UNSIGNED_CHAR ) r=fmrib_main2<char>(argc, argv, src);
    else if ( datatype1==DT_SIGNED_SHORT ) r=fmrib_main2<short>(argc, argv, src);
    else if ( datatype1==DT_SIGNED_INT ) r=fmrib_main2<int>(argc, argv, src);
    else if ( datatype1==DT_FLOAT )  r=fmrib_main2<float>(argc, argv, src);
    else if ( datatype1==DT_DOUBLE ) r=fmrib_main2<double>(argc, argv, src);
	else r=-1;
  }

  freeparser(argc, argv);
  return r;
}

extern "C" __declspec(dllexport) int _stdcall SetaOrigem(char *Volume, char *Saida, short dx, short dy, short dz)
{
  FSLIO *fslio, *dest;
  if (FslFileExists(Volume)) 
  {
    fslio = FslOpen(Volume,"rb");
    dest = FslOpen(Saida, "wb");
	if (fslio != NULL)
	{
        short origin[5]; 
        void *buffer = NULL;
        short x, y, z, v, t;
		int nsrcbytes;

		FslCloneHeader(dest, fslio);
        FslGetDim(fslio, &x, &y, &z, &v);
        nsrcbytes = x * y * z * v * (FslGetDataType(fslio, &t) / 8);

        if( (buffer = calloc(nsrcbytes,1)) == NULL ) {
            return EXIT_FAILURE;
        }
        FslReadVolumes(fslio, buffer, v);
		FslGetAnalyzeOrigin(fslio, origin);

		origin[0] = dx;
		origin[1] = dy;
		origin[2] = dz;

        FslSetAnalyzeSform(dest, origin, fslio->niftiptr->pixdim[1],
                       fslio->niftiptr->pixdim[2], fslio->niftiptr->pixdim[3]);

		FslWriteAllVolumes(dest, buffer);
		FslClose(dest);
		FslFree(dest);

		FslClose(fslio);
		FslFree(fslio);
		free(buffer);
		return 1;
	}
  }
  return 0;
}

extern "C" __declspec(dllexport) int _stdcall PegaOrigem(char *Volume, short *dx, short *dy, short *dz)
{
  FSLIO* fslio;
  if (FslFileExists(Volume)) 
  {
    fslio = FslOpen(Volume,"rb");
	if (fslio != NULL)
	{
        short origin[5]; 
		FslGetAnalyzeOrigin(fslio, origin);

		*dx = origin[0];
		*dy = origin[1];
		*dz = origin[2];

		FslClose(fslio);
		FslFree(fslio);
		return 1;
	}
  }
  return 0;
}

  template <class T>
  int resample(char *innam, float dx, float dy, float dz, int nn, char *Saida)
  {
	string inname=innam;

    volume<T> aniso;
    read_orig_volume(aniso,inname);

    extrapolation oldex = aniso.getextrapolationmethod();
    if ((oldex==boundsassert) || (oldex==boundsexception)) 
      { aniso.setextrapolationmethod(constpad); }

	float stepx, stepy, stepz;
	if (nn==1) aniso.setinterpolationmethod(nearestneighbour);
	else aniso.setinterpolationmethod(trilinear);

	stepx = dx / aniso.xdim();
    stepy = dy / aniso.ydim();
    stepz = dz / aniso.zdim();
    int sx, sy, sz;
    sz = (int) Max(1.0, ( ((float) (aniso.maxz() - aniso.minz() + 1.0)) / stepz));
    sy = (int) Max(1.0, ( ((float) (aniso.maxy() - aniso.miny() + 1.0)) / stepy));
    sx = (int) Max(1.0, ( ((float) (aniso.maxx() - aniso.minx() + 1.0)) / stepx));
    volume<T> iso(sx,sy,sz);
    float fx, fy, fz;
    int x, y, z;
    for (fz=0.0, z=0; z<sz; z++, fz+=stepz) {
      for (fy=0.0, y=0; y<sy; y++, fy+=stepy) {
	for (fx=0.0, x=0; x<sx; x++, fx+=stepx) {
	  iso(x,y,z) = (T)aniso.interpolate(fx,fy,fz);
	}
      }
    }
    iso.copyproperties(aniso);
    iso.setdims(dx,dy,dz);
    // transform the sform and qform matrix appropriately (if set)
    Matrix iso2aniso(4,4);
    iso2aniso = 0.0;
    iso2aniso(1,1)=stepx;
    iso2aniso(2,2)=stepy;
    iso2aniso(3,3)=stepz;
    iso2aniso(4,4)=1.0;
    if (aniso.sform_code()!=NIFTI_XFORM_UNKNOWN) {
      iso.set_sform(aniso.sform_code(), aniso.sform_mat() * iso2aniso);
    }
    if (aniso.qform_code()!=NIFTI_XFORM_UNKNOWN) {
      iso.set_qform(aniso.qform_code(), aniso.qform_mat() * iso2aniso);
    }
	string filename=Saida;
	save_orig_volume(iso, filename);
    aniso.setextrapolationmethod(oldex);
	return 0;
  }


  template <class T>
  int CentralizaVolume(char *Anatomico, char *Saida)
  {
	string inname=Anatomico;
    volume<T> aniso;
    read_orig_volume(aniso,inname);
	int dimmax = Max(aniso.zsize(), Max(aniso.ysize(), aniso.xsize()));
	volume<T> novovolume(dimmax, dimmax, dimmax);
	pad(aniso, novovolume);

	string filename=Saida;
	save_orig_volume(novovolume, filename);
	return 1;
  }

/*  
  template <class T>
  int igualacogs(char *Anatomico, char *Base, char *Saida)
  {
	string inname=Anatomico;
	string basename=Base;
    volume<T> original, base;

	int ox, oy, oz;
	int ax, ay, az;
    read_orig_volume(original,inname);
    read_orig_volume(base,basename);
	Matrix basemat(4,4), origmat(4,4);

    if (base.sform_code()!=NIFTI_XFORM_UNKNOWN) {
	    basemat = base.sform_mat();
    }
    if (base.qform_code()!=NIFTI_XFORM_UNKNOWN) {
	    basemat = base.qform_mat();
    }

    if (original.sform_code()!=NIFTI_XFORM_UNKNOWN) {
	    origmat = original.sform_mat();
    }
    if (original.qform_code()!=NIFTI_XFORM_UNKNOWN) {
	    origmat = original.qform_mat();
    }
	ox = ((float)origmat(1,4)/(float)origmat(1,1));
	oy = ((float)origmat(2,4)/(float)origmat(2,2));
	oz = ((float)origmat(3,4)/(float)origmat(3,3));
 
	ax = ((float)basemat(1,4)/(float)basemat(1,1));
	ay = ((float)basemat(2,4)/(float)basemat(2,2));
	az = ((float)basemat(3,4)/(float)basemat(3,3));

	int offsetx = ax-ox;
	int offsety = ay-oy;
	int offsetz = az-oz;

	volume<T> novovolume(original);
	novovolume.copyproperties(original);
	novovolume=0;
    for (int z=novovolume.minz(); z<=novovolume.maxz(); z++) {
	   for (int y=novovolume.miny(); y<=novovolume.maxy(); y++) {
	      for (int x=novovolume.minx(); x<=novovolume.maxx(); x++) {
	           novovolume(x,y,z) = original(x+offsetx,y+offsety,z+offsetz);
		  }
	   }
	}

    if (original.sform_code()!=NIFTI_XFORM_UNKNOWN) {
	    origmat = original.sform_mat();
		origmat(1,4)=(float)ax/(float)origmat(1,1);
		origmat(2,4)=(float)ay/(float)origmat(2,2);
		origmat(3,4)=(float)az/(float)origmat(3,3);
	    novovolume.set_sform(novovolume.sform_code(), origmat);
    }
    if (original.qform_code()!=NIFTI_XFORM_UNKNOWN) {
	    origmat = original.qform_mat();
		origmat(1,4)=(float)ax/(float)origmat(1,1);
		origmat(2,4)=(float)ay/(float)origmat(2,2);
		origmat(3,4)=(float)az/(float)origmat(3,3);
	    novovolume.set_qform(novovolume.qform_code(), origmat);
    }
	string filename=Saida;
	save_orig_volume(novovolume, filename);
	return 1;
  }

  extern "C" __declspec(dllexport) int _stdcall IgualaCOGs(char *Anatomico, char *Base, char *sAnatomico)
{
  FSLIO* IP1;
  IP1 = FslOpen(Anatomico,"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << Anatomico << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return igualacogs<char>(Anatomico, Base, sAnatomico);
     else if ( datatype==DT_SIGNED_SHORT ) return igualacogs<short>(Anatomico, Base, sAnatomico);
     else if ( datatype==DT_SIGNED_INT ) return igualacogs<int>(Anatomico, Base, sAnatomico);
     else if ( datatype==DT_FLOAT )  return igualacogs<float>(Anatomico, Base, sAnatomico);
     else if ( datatype==DT_DOUBLE ) return igualacogs<double>(Anatomico, Base, sAnatomico);
  }
  return -1;
   
}

  template <class T>
  int DiminuiVolume(char *Anatomico, char *Saida, int xdim, int ydim, int zdim)
  {
	string inname=Anatomico;
    volume<T> original;
    read_orig_volume(original,inname);
	int tamx =  original.maxx() - original.minx();
	int tamy =  original.maxy() - original.miny();
	int tamz =  original.maxz() - original.minz();

	int offsetx = tamx-xdim;
	int offsety = tamy-ydim;
	int offsetz = tamz-zdim;

	volume<T> novovolume(xdim, ydim, zdim);
	novovolume.copyproperties(original);
    for (int z=novovolume.minz(); z<=novovolume.maxz(); z++) {
	   for (int y=novovolume.miny(); y<=novovolume.maxy(); y++) {
	      for (int x=novovolume.minx(); x<=novovolume.maxx(); x++) {
	           novovolume(x,y,z) = original(x+offsetx,y+offsety,z+offsetz);
		  }
	   }
	}

	Matrix pad2vol(4,4);
    pad2vol = IdentityMatrix(4);
    pad2vol(1,4) = offsetx;
    pad2vol(2,4) = offsety;
    pad2vol(3,4) = offsetz;
    if (novovolume.sform_code()!=NIFTI_XFORM_UNKNOWN) {
	    novovolume.set_sform(novovolume.sform_code(), novovolume.sform_mat() * pad2vol);
    }
    if (novovolume.qform_code()!=NIFTI_XFORM_UNKNOWN) {
	    novovolume.set_qform(novovolume.qform_code(),novovolume.qform_mat() * pad2vol);
    }

	string filename=Saida;
	save_orig_volume(novovolume, filename);
	return 1;
  }

extern "C" __declspec(dllexport) int _stdcall ShrinkVolume(char *Anatomico, char *sAnatomico, int newxdim, int newydim, int newzdim)
{
  FSLIO* IP1;
  IP1 = FslOpen(Anatomico,"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << Anatomico << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return DiminuiVolume<char>(Anatomico, sAnatomico, newxdim, newydim, newzdim);
     else if ( datatype==DT_SIGNED_SHORT ) return DiminuiVolume<short>(Anatomico, sAnatomico, newxdim, newydim, newzdim);
     else if ( datatype==DT_SIGNED_INT ) return DiminuiVolume<int>(Anatomico, sAnatomico, newxdim, newydim, newzdim);
     else if ( datatype==DT_FLOAT )  return DiminuiVolume<float>(Anatomico, sAnatomico, newxdim, newydim, newzdim);
     else if ( datatype==DT_DOUBLE ) return DiminuiVolume<double>(Anatomico, sAnatomico, newxdim, newydim, newzdim);
  }
  return -1;
   
}
*/
extern "C" __declspec(dllexport) int _stdcall Realinha(char *Anatomico, char *sAnatomico)
{
  FSLIO* IP1;
  IP1 = FslOpen(Anatomico,"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << Anatomico << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
     if ( datatype==DT_UNSIGNED_CHAR ) return CentralizaVolume<char>(Anatomico, sAnatomico);
     else if ( datatype==DT_SIGNED_SHORT ) return CentralizaVolume<short>(Anatomico, sAnatomico);
     else if ( datatype==DT_SIGNED_INT ) return CentralizaVolume<int>(Anatomico, sAnatomico);
     else if ( datatype==DT_FLOAT )  return CentralizaVolume<float>(Anatomico, sAnatomico);
     else if ( datatype==DT_DOUBLE ) return CentralizaVolume<double>(Anatomico, sAnatomico);
  }
  return -1;
   
}

extern "C" __declspec(dllexport) int _stdcall AnatomicoResampleado(char *Anatomico, char *sAnatomico, float dx, float dy, float dz, float TR, int nn)
{
/*
  FSLIO* fslio;
  char CmdLn[1000];
  if (FslFileExists(Anatomico)) 
  {
    fslio = FslOpen(Anatomico,"rb");
	if (fslio != NULL)
	{
        short origin[5]; 
		short x, y, z, v, ox, oy, oz;
		short nx, ny, nz;
		FslGetAnalyzeOrigin(fslio, origin);
        FslGetDim(fslio,&x,&y,&z,&v);
		nx = (short) floor((x / dx) + 0.5);
		ny = (short) floor((y / dy) + 0.5);
		nz = (short) floor((z / dz) + 0.5);

		ox = (short) floor((origin[0] / dx) + 0.5);
		oy = (short) floor((origin[1] / dy) + 0.5);
		oz = (short) floor((origin[2] / dz) + 0.5);
        sprintf(CmdLn, "fslcreatehd %d %d %d 1 %f %f %f %f %d %d %d 2 %s", nx, ny, nz, dx, dy, dz, TR, ox, oy, oz, sAnatomico);
		FslClose(fslio);
		FslFree(fslio);
		return fslcreatehd(CmdLn);
	}
  }
*/
  FSLIO* IP1;
  IP1 = FslOpen(Anatomico,"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << Anatomico << " for reading!\n";
  }
  else
  {
     short tipo;
     FslGetDataType(IP1,&tipo);
     FslClose(IP1);
     FslFree(IP1);

     short datatype=NEWIMAGE::closestTemplatedType(tipo);
	 if (nn ==0) return resample<float>(Anatomico, dx, dy, dz, nn, sAnatomico);
	 else
	 {
        if ( datatype==DT_UNSIGNED_CHAR ) return resample<char>(Anatomico, dx, dy, dz, nn, sAnatomico);
        else if ( datatype==DT_SIGNED_SHORT ) return resample<short>(Anatomico, dx, dy, dz, nn, sAnatomico);
        else if ( datatype==DT_SIGNED_INT ) return resample<int>(Anatomico, dx, dy, dz, nn, sAnatomico);
        else if ( datatype==DT_FLOAT )  return resample<float>(Anatomico, dx, dy, dz, nn, sAnatomico);
        else if ( datatype==DT_DOUBLE ) return resample<double>(Anatomico, dx, dy, dz, nn, sAnatomico);
	 }
  }
  return -1;
   
}

extern "C" __declspec(dllexport) void _stdcall CalculaMedias(char *arquivobase, char*arqintervalos)
{
	Intervalo *intervalos;
	char **condicoes;
	int tamintervalos, tamcondicoes;
	char DirSaida[255], media[255], temp[255];

	extractfilepath(arquivobase, DirSaida);
	if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");

	sprintf(media, "%s%s", DirSaida, "media");
	sprintf(temp, "%s%s", DirSaida, "temp.nii");

    tamintervalos=LeIntervalo(arqintervalos, &intervalos);
    tamcondicoes=pegalistacondicoes(intervalos, tamintervalos, &condicoes);
	for (int t=0; t<tamintervalos; t++)
	{
		char cmd[500], arqmedia[500];
		sprintf(arqmedia, "%s%d.nii", media, t+1);
		if (!fileexists(arqmedia))
		{
			sprintf(cmd, "fslroi %s %s %d %d", arquivobase, temp, intervalos[t].inicio-1, intervalos[t].fim-intervalos[t].inicio+1);
			fslroi(cmd);
			sprintf(cmd, "fslmaths %s -Tmean %s%d.nii", temp, media, t+1);
			fslmaths(cmd);
		}
	}
	free(intervalos);
	desalocacondicoes(condicoes, tamcondicoes);
	remove(temp);
}

extern "C" __declspec(dllexport) void _stdcall NormalizacaoZ(char *arquivobase, char *mascara, char *saida)
{
	FSLIO *arq;
	char cmd[2000];
	char meanfile[500], stddevfile[500];
	char DirSaida[500];
	int total;

	extractfilepath(arquivobase, DirSaida);
	if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");

	arq = fslioopen(arquivobase);
	total=arq->niftiptr->nt;
	fslioclose(arq);

	sprintf(meanfile, "%s%s", DirSaida, "mean.nii");
	sprintf(stddevfile, "%s%s", DirSaida, "stddev.nii");
//	if (!fileexists(meanfile)) 
	{
		sprintf(cmd, "fslmaths %s -Tmean %s", arquivobase, meanfile);
		fslmaths(cmd);
    }

//	if (!fileexists(stddevfile)) 
	{
		sprintf(cmd, "fslmaths %s -Tstd %s", arquivobase, stddevfile);
		fslmaths(cmd);
	}

	sprintf(cmd, "fslmaths %s -sub %s -div %s %s", arquivobase, meanfile, stddevfile, saida);
	fslmaths(cmd);
		   
    if (mascara != NULL)
    {
	   sprintf(cmd, "fslmaths %s -mas %s %s", saida, mascara, saida);
	   fslmaths(cmd);
	}
	remove(meanfile);
	remove(stddevfile);
}

extern "C" __declspec(dllexport) void _stdcall PlotaProjecoes(char *arquivo, char *xlabel, char *ylabel, char *saida)
{
  Matrix in, fake(100,1);
  fake = 0;
  string Arquivo = arquivo;
  string Saida = saida;
  in=read_ascii_matrix(Arquivo);

  miscplot newplot;

  newplot.set_minmaxscale(1.001);
//  for (int i = 0; i < (int)labels.value().size(); i++) newplot.add_label(labels.value().at(i));
  newplot.set_xysize(in.Minimum(),in.Maximum());

  newplot.set_yrange(1,5);
  newplot.add_xlabel(xlabel);
  newplot.add_ylabel(ylabel); 
  newplot.setscatter(in, 2);
  newplot.timeseries(fake, Saida, "Projections");
}

extern "C" __declspec(dllexport) void _stdcall NormalizacaoZ2(char *arquivobase, char *mascara, char *saida)
{
	FSLIO *arq;
	char cmd[2000];
	char meanfile[500], stddevfile[500];
	char DirSaida[500];
	int total;

	extractfilepath(arquivobase, DirSaida);
	if (strlen(DirSaida) > 0) strcat(DirSaida, "\\");

	arq = fslioopen(arquivobase);
	total=arq->niftiptr->nt;
	fslioclose(arq);

	sprintf(meanfile, "%s%s", DirSaida, "mean.nii");
	sprintf(stddevfile, "%s%s", DirSaida, "stddev.nii");

	sprintf(cmd, "fslmaths %s -sub %s -div %s %s", arquivobase, meanfile, stddevfile, saida);
	fslmaths(cmd);
		   
	sprintf(cmd, "fslmaths %s -mas %s %s", saida, mascara, saida);
	fslmaths(cmd);
}

extern "C" __declspec(dllexport) void _stdcall CriaVolume4D(volume<float> *ref, int tamanho, volume4D<float> *&saida)
{
	saida = new volume4D<float>;
	saida->reinitialize(ref->xsize(), ref->ysize(), ref->zsize(), tamanho);
	saida->copyproperties(*ref);
}

extern "C" __declspec(dllexport) void _stdcall Subtracao4D(volume4D<float> *vol4D, int ini, int fim, volume<float> *aretirar)
{
	if (vol4D->tsize() ==1) (*vol4D)[0]-=(*aretirar);
	else
	{
		for (int t=ini; t<=fim; t++)
			(*vol4D)[t-ini]-=(*aretirar);
	}
}

extern "C" __declspec(dllexport) void* _stdcall PonteiroVolume4D(volume4D<float> *vol4D, int indice)
{
	return &(*vol4D)[indice-1](0, 0, 0);
}

extern "C" __declspec(dllexport) void* _stdcall PonteiroVolume(volume<float> *vol)
{
	return &(*vol)(0, 0, 0);
}

extern "C" __declspec(dllexport) void _stdcall AdicionaVolume4D(volume4D<float> *vol4D, int indice, volume4D<float> *aadicionar)
{
	(*vol4D)[indice-1].copydata((*aadicionar)[0]);
}

extern "C" __declspec(dllexport) void _stdcall DuplicaVolume4DIndice(volume4D<float> *vol4D, int indice, volume4D<float> *&duplicado)
{
	duplicado = new volume4D<float>;
	duplicado->reinitialize(vol4D->xsize(), vol4D->ysize(), vol4D->zsize(), 1, 0);
	duplicado->copyproperties(*vol4D);
	(*duplicado)[0].copydata((*vol4D)[indice-1]);
}

extern "C" __declspec(dllexport) void _stdcall DuplicaVolume(volume<float> *vol, volume<float> *&duplicado)
{
	duplicado = new volume<float>;
	duplicado->reinitialize(*vol);
	duplicado->copyproperties(*vol);
	duplicado->copydata(*vol);
}

extern "C" __declspec(dllexport) void _stdcall DuplicaVolume4D(volume4D<float> *vol4D, volume4D<float> *&duplicado)
{
	duplicado = new volume4D<float>;
	duplicado->reinitialize(*vol4D);
	duplicado->copyproperties(*vol4D);
	for(int t=0; t<vol4D->tsize(); t++) (*duplicado)[t].copydata((*vol4D)[t]);
}

extern "C" __declspec(dllexport) void _stdcall DuplicaVolume3D(volume<float> *vol, volume4D<float> *&duplicado)
{
	duplicado = new volume4D<float>;
	duplicado->reinitialize(vol->xsize(), vol->ysize(), vol->zsize(), 1);
	duplicado->copyproperties(*vol);
	(*duplicado)[0].copydata(*vol);
}

extern "C" __declspec(dllexport) void _stdcall FiltroGaussiano(volume4D<float> *vol, float FWHM)
{
   float xdim=vol->xdim();
   float ydim=vol->ydim();
   float zdim=vol->zdim();
   bool separable(true);
   volume<float> kernel(box_kernel(3,3,3));
   kernel=gaussian_kernel3D(FWHM,xdim,ydim,zdim);
   *vol=generic_convolve(*vol,kernel,separable,true);
}

extern "C" __declspec(dllexport) void _stdcall SetaTipoSaida(int tipo)
{
	FslSetOverrideOutputType(tipo);
}

extern "C" __declspec(dllexport) void _stdcall ParametrosMatriz(char *matriz, char *volref, char *saida)
{
   float rmax = 80.0;
   ColumnVector param_vec(12);
   ColumnVector center(3);
   Matrix Matriz = read_ascii_matrix(string(matriz));

   volume<float>refvol;
   read_volume(refvol, string(volref));

   ofstream outfile, rmsabsfile;
   string filename = string(saida) + ".par";
   string rms_abs_filename = string(saida) + "_abs.rms";

   param_vec = 0;   
   center(1) = 0.5*(refvol.xsize() - 1.0)*refvol.xdim();
   center(2) = 0.5*(refvol.ysize() - 1.0)*refvol.ydim();
   center(3) = 0.5*(refvol.zsize() - 1.0)*refvol.zdim();
   
   rmsabsfile. open(rms_abs_filename.c_str());
   rmsabsfile << rms_deviation(IdentityMatrix(4), Matriz, center, rmax) << endl;
 
   outfile. open(filename.c_str());
   decompose_aff(param_vec, Matriz, refvol.cog("scaled_mm"), rotmat2euler);
   outfile << param_vec(1) << "  " << param_vec(2) << "  " 
		 << param_vec(3) << "  " << param_vec(4) << "  " 
		 << param_vec(5) << "  " << param_vec(6) << "  " << endl;
}