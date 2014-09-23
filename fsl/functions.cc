#include "miscmaths/miscprob.h"
#include "miscmaths/t2z.h"
#include "libprob.h"

#include <cstring>

using namespace MISCMATHS;
using namespace NEWIMAGE;

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

char *find_line(char *filename, char *key, char *fl)
{
  FILE *fd;
  char *return_ptr=NULL, tmp_fl[10000];
  int  j;

  fd=fopen(filename,"rb");

  while ( fgets(tmp_fl, 1000, fd) != NULL )
    if (strncmp(tmp_fl,"set ",4)==0)
      for(j=4; tmp_fl[j]!=0; j++)
	if (strncmp(key,tmp_fl+j,strlen(key))==0)
	  {
	    strcpy(fl,tmp_fl+j);
	    return_ptr = fl+1+strlen(key);
	    if (fl[strlen(fl)-1]==10) fl[strlen(fl)-1]=0;
	    if (fl[strlen(fl)-1]==13) fl[strlen(fl)-1]=0;
	    if (fl[strlen(fl)-1]=='"') fl[strlen(fl)-1]=0;
	    if (fl[1+strlen(key)]=='"') return_ptr++;
	  }

  fclose(fd);

  if (return_ptr==NULL)
    {
      printf("Can't find key %s\n",key);
      exit(1);
    }
  else
    return return_ptr;
}


extern "C" __declspec(dllexport) int _stdcall CalculaAtivacao(int Inicio, int ini, int fim, int tamb, char *prefixo, char *Basal, char *Mascara, char *Saida)
{
   char tmpname[5000];

   volume<float> basal;
   if (Basal != NULL) 
   {
      string bName = Basal;
      read_volume(basal, bName);
   }
   
   volume<float> auxvol;
   sprintf(tmpname, "%s%.4d.nii", prefixo, ini);
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
	      sprintf(tmpname, "%s%.4d.nii", prefixo, j);
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

extern "C" __declspec(dllexport) int _stdcall convolui(char *matrixfile, char *volnam, float *rx, float *ry, float *rz, float *tx, float *ty, float *tz)
{
	float sigma1=2.449, delay1=6, // first gamma
	sigma2=4, delay2=16,    // second gamma
	ratio=6;                // hrf = gammapdf1 - gammapdf2/ratio;

	int fw = (int)((delay2 + sigma2*5)*mult);
	maxconvwin=MAX(fw,maxconvwin);

	ColumnVector cX = mygammapdf(fw,mult,delay1,sigma1) - mygammapdf(fw,mult,delay2,sigma2)/ratio;

	ColumnVector oX=do_convolve(orig_X.Column(orig_ev+1),cX,convolve_phase,1);
	do_resample(oX, real_X, real_ev, trmult, negpts[orig_ev]);
}

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
	 FslClose(fslio);

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
  
     invol.swapdimensions(newx,newy,newz);

     if (outname!="") {
        retval = save_orig_volume4D(invol,outname);
     }
  }
  return retval;
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
     free(IP1);

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


  if (affmat.Determinant()<0.0) {
    cout << "WARNING:: Flipping Left/Right orientation (as det < 0)" << endl;
  }
  
  invol.swapdimensions(newx,newy,newz);

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

extern "C" __declspec(dllexport) int _stdcall fslswapdim_rt(char *CmdLn, FSLIO *src)
{
  int r;
  int argc;
  char **argv;
  
  parser(CmdLn, argc, argv);

  string progname=argv[0];
  if (argc<5) { 
    print_usage(progname);
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
