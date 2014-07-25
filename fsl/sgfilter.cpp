#include <math.h> 
#include <stdlib.h> 
#include <stdarg.h> 
#include <stdio.h> 
#include <string.h> 
#include <ctype.h> 
#include <time.h> 
#include "sgfilter.h"
#include <vector>

using namespace std;

extern char*optarg;
char*progname;

int* ivector(long nl,long nh)
{
   int*v;
   v = (int*)malloc((size_t)((nh-nl+2)*sizeof(int)));
   if(!v)
   {
      printf("Error: Allocation failure.");
      exit(1);
   }
   return v-nl+1;
}

double* dvector(long nl,long nh)
{
   double*v;
   long k;
   v = (double*)malloc((size_t)((nh-nl+2)*sizeof(double)));
   if(!v)
   {
      printf("Error: Allocation failure.");
      exit(1);
   }
   for(k=nl; k<=nh; k++) v[k]= 0.0;
   return v-nl+1;
}

double** dmatrix(long nrl,long nrh,long ncl,long nch)
{
   long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double**m;
   m=(double**)malloc((size_t)((nrow+1)*sizeof(double*)));
   if(!m)
   {
      printf("Allocation failure 1 occurred.");
      exit(1);
   }
   m+= 1;
   m-= nrl;
   m[nrl]= (double*)malloc((size_t)((nrow*ncol+1)*sizeof(double)));
   if(!m[nrl])
   {
      printf("Allocation failure 2 occurred.");
      exit(1);
   }
   m[nrl]+= 1;
   m[nrl]-= ncl;
   for(i= nrl+1;i<=nrh;i++)m[i]= m[i-1]+ncol;
   return m;
}

void free_ivector(int*v,long nl,long nh)
{
   free((char*)(v+nl-1));
}

void free_dvector(double*v,long nl,long nh)
{
   free((char*)(v+nl-1));
}

void free_dmatrix(double**m,long nrl,long nrh,long ncl,long nch)
{
   free((char*)(m[nrl]+ncl-1));
   free((char*)(m+nrl-1));
}

void lubksb(double**a,int n,int*indx,double b[])
{
   int i,ii= 0,ip,j;
   double sum;

   for(i= 1;i<=n;i++)
   {
      ip= indx[i];
      sum= b[ip];
      b[ip]= b[i];
      if(ii)
         for(j= ii;j<=i-1;j++)sum-= a[i][j]*b[j];
      else if(sum)ii= i;
      b[i]= sum;
   }
   for(i= n;i>=1;i--)
   {
      sum= b[i];
      for(j= i+1;j<=n;j++)sum-= a[i][j]*b[j];
      b[i]= sum/a[i][i];
   }
}

void ludcmp(double**a,int n,int*indx,double*d)
{
   int i,imax= 0,j,k;
   double big,dum,sum,temp;
   double*vv;

   vv= dvector(1,n);
   *d= 1.0;
   for(i= 1;i<=n;i++)
   {
      big= 0.0;
      for(j= 1;j<=n;j++)
      if((temp= fabs(a[i][j]))> big) big = temp;
      if(big==0.0)
	  {
         printf("Error: Singular matrix found in routine ludcmp()");
         exit(1);
      }
      vv[i]= 1.0/big;
   }
   for(j= 1;j<=n;j++)
   {
      for(i= 1;i<j;i++)
	  {
         sum= a[i][j];
         for(k= 1;k<i;k++)sum-= a[i][k]*a[k][j];
         a[i][j]= sum;
      }
      big= 0.0;
      for(i= j;i<=n;i++)
	  {
         sum= a[i][j];
         for(k= 1;k<j;k++)
            sum-= a[i][k]*a[k][j];
         a[i][j]= sum;
         if((dum= vv[i]*fabs(sum))>=big)
		 {
            big = dum;
            imax = i;
         }
      }
      if(j!=imax)
	  {
         for(k= 1;k<=n;k++)
		 {
            dum= a[imax][k];
            a[imax][k]= a[j][k];
            a[j][k]= dum;
         }
         *d= -(*d);
         vv[imax]= vv[j];
      }
      indx[j]= imax;
      if(a[j][j]==0.0)a[j][j]= EPSILON;
      if(j!=n)
	  {
         dum= 1.0/(a[j][j]);
         for(i= j+1;i<=n;i++)a[i][j]*= dum;
      }
   }
   free_dvector(vv,1,n);
}

#define SWAP(a,b) tempr= (a);(a)= (b);(b)= tempr
void four1(double data[],unsigned long nn,int isign)
{
   unsigned long n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;
   double tempr,tempi;

   n= nn<<1;
   j= 1;
   for(i= 1;i<n;i+= 2)
   {
      if(j> i)
	  {
         SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
      }
      m= n>>1;
      while(m>=2&&j> m)
	  {
         j-= m;
         m>>= 1;
      }
      j+= m;
   }
   mmax= 2;
   while(n> mmax)
   {
      istep= mmax<<1;
      theta= isign*(6.28318530717959/mmax);
      wtemp= sin(0.5*theta);
      wpr= -2.0*wtemp*wtemp;
      wpi= sin(theta);
      wr= 1.0;
      wi= 0.0;
      for(m= 1;m<mmax;m+= 2)
	  {
         for(i= m;i<=n;i+= istep)
		 {
            j= i+mmax;
            tempr= wr*data[j]-wi*data[j+1];
            tempi= wr*data[j+1]+wi*data[j];
            data[j]= data[i]-tempr;
            data[j+1]= data[i+1]-tempi;
            data[i]+= tempr;
            data[i+1]+= tempi;
         }
         wr= (wtemp= wr)*wpr-wi*wpi+wr;
         wi= wi*wpr+wtemp*wpi+wi;
      }
      mmax= istep;
   }
}
#undef SWAP

void twofft(double data1[], double data2[], double fft1[], double fft2[], unsigned long n)
{
   unsigned long nn3,nn2,jj,j;
   double rep,rem,aip,aim;

   nn3= 1+(nn2= 2+n+n);
   for(j= 1,jj= 2;j<=n;j++,jj+= 2)
   {
      fft1[jj-1]= data1[j];
      fft1[jj]= data2[j];
   }
   four1(fft1,n,1);
   fft2[1]= fft1[2];
   fft1[2]= fft2[2]= 0.0;
   for(j= 3;j<=n+1;j+= 2)
   {
      rep= 0.5*(fft1[j]+fft1[nn2-j]);
      rem= 0.5*(fft1[j]-fft1[nn2-j]);
      aip= 0.5*(fft1[j+1]+fft1[nn3-j]);
      aim= 0.5*(fft1[j+1]-fft1[nn3-j]);
      fft1[j]= rep;
      fft1[j+1]= aim;
      fft1[nn2-j]= rep;
      fft1[nn3-j]= -aim;
      fft2[j]= aip;
      fft2[j+1]= -rem;
      fft2[nn2-j]= aip;
      fft2[nn3-j]= rem;
   }
}

void realft(double data[],unsigned long n,int isign)
{
   unsigned long i,i1,i2,i3,i4,np3;
   double c1= 0.5,c2,h1r,h1i,h2r,h2i;
   double wr,wi,wpr,wpi,wtemp,theta;

   theta= 3.141592653589793/(double)(n>>1);
   if(isign==1)
   {
      c2= -0.5;
      four1(data,n>>1,1);
   }
   else
   {
      c2= 0.5;
      theta= -theta;
   }
   wtemp= sin(0.5*theta);
   wpr= -2.0*wtemp*wtemp;
   wpi= sin(theta);
   wr= 1.0+wpr;
   wi= wpi;
   np3= n+3;
   for(i= 2;i<=(n>>2);i++)
   {
      i4= 1+(i3= np3-(i2= 1+(i1= i+i-1)));
      h1r= c1*(data[i1]+data[i3]);
      h1i= c1*(data[i2]-data[i4]);
      h2r= -c2*(data[i2]+data[i4]);
      h2i= c2*(data[i1]-data[i3]);
      data[i1]= h1r+wr*h2r-wi*h2i;
      data[i2]= h1i+wr*h2i+wi*h2r;
      data[i3]= h1r-wr*h2r+wi*h2i;
      data[i4]= -h1i+wr*h2i+wi*h2r;
      wr= (wtemp= wr)*wpr-wi*wpi+wr;
      wi= wi*wpr+wtemp*wpi+wi;
   }
   if(isign==1)
   {
      data[1]= (h1r= data[1])+data[2];
      data[2]= h1r-data[2];
   }
   else
   {
      data[1]= c1*((h1r= data[1])+data[2]);
      data[2]= c1*(h1r-data[2]);
      four1(data,n>>1,-1);
   }
}

char convlv(double data[],unsigned long n,double respns[],unsigned long m,
int isign,double ans[])
{
   unsigned long i,no2;
   double dum,mag2,*fft;

   fft= dvector(1,n<<1);

   for(i= 1;i<=(m-1)/2;i++)
      respns[n+1-i]= respns[m+1-i];

   for(i= (m+3)/2;i<=n-(m-1)/2;i++)
      respns[i]= 0.0;

   twofft(data,respns,fft,ans,n);
   no2= n>>1;
   for(i= 2;i<=n+2;i+= 2)
   {
      if(isign==1)
	  {
         ans[i-1]= (fft[i-1]*(dum= ans[i-1])-fft[i]*ans[i])/no2;
         ans[i]= (fft[i]*dum+fft[i-1]*ans[i])/no2;
      }
	  else 
		  if(isign==-1)
		  {
             if((mag2= ans[i-1]*ans[i-1]+ans[i]*ans[i])==0.0)
			 {
                printf("Attempt of deconvolving at zero response in convlv().");
                return(1);
             }
             ans[i-1]= (fft[i-1]*(dum= ans[i-1])+fft[i]*ans[i])/mag2/no2;
             ans[i]= (fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
          }
		  else
		  {
             printf("No meaning for isign in convlv().");
             return(1);
          }
   }
   ans[2]= ans[n+1];
   realft(ans,n,-1);
   free_dvector(fft,1,n<<1);
   return(0);
}

char sgcoeff(double c[],int np,int nl,int nr,int ld,int m)
{
   int imj,ipj,j,k,kk,mm,*indx;
   double d,fac,sum,**a,*b;

   if(np<nl+nr+1||nl<0||nr<0||ld> m||nl+nr<m)
   {
      printf("Inconsistent arguments detected in routine sgcoeff.");
      return(1);
   }
   indx= ivector(1,m+1);
   a= dmatrix(1,m+1,1,m+1);
   b= dvector(1,m+1);
   for(ipj= 0;ipj<=(m<<1);ipj++)
   {
      sum= (ipj?0.0:1.0);

      for (k= 1;k<=nr;k++)
		  sum+= pow((double)k,(double)ipj);

      for (k= 1;k<=nl;k++)
		  sum+= pow((double)-k,(double)ipj);

      mm= (ipj<2*m-ipj?ipj:2*m-ipj);
      for (imj= -mm;imj<=mm;imj+= 2) 
		  a[1+(ipj+imj)/2][1+(ipj-imj)/2]= sum;
   }
   ludcmp(a,m+1,indx,&d);
   for(j= 1;j<=m+1;j++)
	   b[j]= 0.0;
   b[ld+1]= 1.0;
   lubksb(a,m+1,indx,b);
   for(kk= 1;kk<=np;kk++)
	   c[kk]= 0.0;
   for(k= -nl;k<=nr;k++)
   {
      sum= b[1];
      fac= 1.0;
      for(mm= 1;mm<=m;mm++)
		  sum+= b[mm+1]*(fac*= k);
      kk= ((np-k)%np)+1;
      c[kk]= sum;
   }
   free_dvector(b,1,m+1);
   free_dmatrix(a,1,m+1,1,m+1);
   free_ivector(indx,1,m+1);
   return(0);
}

char sgfilter(vector<double> &yr, vector<double> &yf, int nl, int nr, int ld, int m)
{
   int np= nl+1+nr;
   double*c;
   char retval;

#if CONVOLVE_WITH_NR_CONVLV 
   c= dvector(1,mm);
   retval= sgcoeff(c,np,nl,nr,ld,m);
   if(retval==0)
   convlv(yr,mm,c,np,1,yf);
   free_dvector(c,1,mm);
#else 
   int j;
   long int k;

   yf.resize(yr.size());
   c= dvector(1,nl+nr+1);
   retval= sgcoeff(c, np, nl, nr, ld, m);
   if(retval==0)
   {
      for(k= 1;k<=nl;k++)
	  {
         for(yf[k-1]= 0.0,j= -nl;j<=nr;j++)
		 {
            if(k+j>=1)
			{
               yf[k-1]+= c[(j>=0?j+1:nr+nl+2+j)] * yr[k+j-1];
            }
         }
      }
      for(k= nl+1;k<=yr.size()-nr;k++)
	  {
         for(yf[k-1]= 0.0,j= -nl;j<=nr;j++)
		 {
            yf[k-1]+= c[(j>=0?j+1:nr+nl+2+j)]*yr[k+j-1];
         }
      }
      for(k=yr.size()-nr+1;k <= yr.size();k++)
	  {
         for(yf[k-1]= 0.0,j= -nl;j<=nr;j++)
		 {
            if(k+j <= yr.size())
			{
               yf[k-1]+= c[(j>=0?j+1:nr+nl+2+j)]*yr[k+j-1];
            }
         }
      }
   }
   free_dvector(c,1,nr+nl+1);
#endif
   return(retval);
}

short pathcharacter(int ch)
{
   return(isalnum(ch)||(ch=='.')||(ch=='/')||(ch=='\\')||(ch=='_')||
    (ch=='-')||(ch=='+'));
}

char*strip_away_path(char filename[])
{
   int j,k= 0;

   while(pathcharacter(filename[k]))
	   k++;

   j= (--k);

   while(isalnum((int)(filename[j])))
	   j--;
   j++;
   return(&filename[j]);
}

void hl(const char*format,...)
{
   va_list args;
   char line[1024];

   va_start(args,format);
   vsprintf(line,format,args);
   va_end(args);
   sprintf(line+strlen(line),"\n");
   fprintf(stdout,"%s",line);
   return;
}

void showsomehelp(void)
{
   hl("Usage: %s [options]",progname);
   hl("Options:");
   hl(" -h, --help");
   hl("    Display this help message and exit clean.");
   hl(" -i, --inputfile <str>");
   hl("    Specifies the file name from which unfiltered data is to be read.");
   hl("    The input file should describe the input as two columns contain-");
   hl("    ing $x$- and $y$-coordinates of the samples.");
   hl(" -o, --outputfile <str>");
   hl("    Specifies the file name to which filtered data is to be written,");
   hl("    again in a two-column format containing $x$- and $y$-coordinates");
   hl("    of the filtered samples. If this option is omitted, the generated");
   hl("    filtered data will instead be written to the console (terminal).");
   hl(" -nl <nl>");
   hl("    Specifies the number of samples nl to use to the 'left' of the");
   hl("    basis sample in the regression window (kernel). The total number");
   hl("    of samples in the window will be nL+nR+1.");
   hl(" -nr <nr>");
   hl("    Specifies the number of samples nr to use to the 'right' of the");
   hl("    basis sample in the regression window (kernel). The total number");
   hl("    of samples in the window will be nL+nR+1.");
   hl(" -m <m>");
   hl("    Specifies the order m of the polynomial to use in the regression");
   hl("    analysis leading to the Savitzky-Golay coefficients. Typical");
   hl("    values are between m=2 and m=6. Beware of too high values, which");
   hl("    easily makes the regression too sensitive, with an oscillatory");
   hl("    result.");
   hl(" -ld <ld>");
   hl("    Specifies the order of the derivative to extract from the ");
   hl("    Savitzky--Golay smoothing algorithm. For regular Savitzky-Golay");
   hl("    smoothing of the input data as such, use ld=0. For the Savitzky-");
   hl("    Golay smoothing and extraction of derivatives, set ld to the");
   hl("    order of the desired derivative and make sure that you correctly");
   hl("    interpret the scaling parameters as described in 'Numerical");
   hl("    Recipes in C', 2nd Edn (Cambridge University Press, New York,");
   hl("    1994).");
   hl(" -v, --verbose");
   hl("    Toggle verbose mode. (Default: Off.)  This option should always");
   hl("    be omitted whenever no  output file has been specified (that is");
   hl("    to say, omit any --verbose or -v option whenever --outputfile or");
   hl("    -o has been omitted), as the verbose logging otherwise will");
   hl("    contaminate the filtered data stream written to the console");
   hl("    (terminal).");
}

long int num_coordinate_pairs(FILE*file)
{
   double tmp;
   int tmpch;
   long int mm= 0;
   fseek(file,0L,SEEK_SET);
   while((tmpch= getc(file))!=EOF)
   {
      ungetc(tmpch,file);
      fscanf(file,"%lf",&tmp);
      //fscanf(file,"%lf",&tmp);
      mm++;
      tmpch= getc(file);
      while((tmpch!=EOF)&&(!isdigit(tmpch)))
		  tmpch= getc(file);

      if(tmpch!=EOF)ungetc(tmpch,file);
   }
   fseek(file,0L,SEEK_SET);
   return(mm);
}

extern "C" __declspec(dllexport) void _stdcall retornaComponenteDetrending(int tamanho, float *timeserie, float *componente, int nl, int nr, int m, int ld)
{
	vector<double> yr, yf;

	yr.resize(tamanho);
	yf.resize(tamanho);

	for (int i = 0; i<tamanho;i++) yr[i] = timeserie[i];

    sgfilter(yr, yf, nl, nr, ld, m);

	for (int i = 0; i<tamanho;i++) componente[i] = yf[i];
}