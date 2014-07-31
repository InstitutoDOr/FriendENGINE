#define VERSION "1.6"
#define COPYRIGHT "Copyright (C) 2006-2011, Fredrik Jonsson"
#define DEFAULT_NL (15)
#define DEFAULT_NR (15)
#define DEFAULT_M (4)
#define DEFAULT_LD (0)
#define EPSILON ((double)(1.0e-20))
#define NCHMAX (256)
#define CONVOLVE_WITH_NR_CONVLV (0)


#if defined(_CYGWIN_SIGNAL_H)||defined(__APPLE__)||defined(__unix__)||defined(__linux)
#define UNIX_LIKE_OS (1)
#endif

int*ivector(long nl,long nh);
double*dvector(long nl,long nh);
double**dmatrix(long nrl,long nrh,long ncl,long nch);
void free_ivector(int*v,long nl,long nh);
void free_dvector(double*v,long nl,long nh);
void free_dmatrix(double**m,long nrl,long nrh,long ncl,long nch);
void lubksb(double**a,int n,int*indx,double b[]);
void ludcmp(double**a,int n,int*indx,double*d);
void four1(double data[],unsigned long nn,int isign);
void twofft(double data1[],double data2[],double fft1[],double fft2[],
unsigned long n);
void realft(double data[],unsigned long n,int isign);
char convlv(double data[],unsigned long n,double respns[],unsigned long m,
int isign,double ans[]);
char sgcoeff(double c[],int np,int nl,int nr,int ld,int m);
char sgfilter(double yr[],double yf[],int mm,int nl,int nr,int ld,int m);
char*strip_away_path(char filename[]);
long int num_coordinate_pairs(FILE*file);
