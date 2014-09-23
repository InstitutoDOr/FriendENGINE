// {{{ Copyright etc. 

/*  feat_model - create FEAT design matrix, contrasts etc.

    Stephen Smith and Matthew Webster, FMRIB Analysis Group

    Copyright (C) 1999-2009 University of Oxford  */
/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 4.0 (c) 2007, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/1112. */

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "featlib.h"

#include "miscmaths/miscprob.h"
#include "miscmaths/t2z.h"
#include "libprob.h"

#include <cstring>
#include "parser.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;

#define DT                0.05 /* temporal sampling for initial model, in seconds */
#define NEGSECS           30.0 /* amount of negative modelling allowed, for custom 3 etc, in seconds */

// }}}
// {{{ z2t 

namespace featmodel {
extern "C" __declspec(dllexport) int _stdcall wpng(char *CmdLn);

float z2t(float z, int dof)
{
  float tmin=0, tmax=1e10, absz=fabs(z), old_diff=1e11;
  
  while (old_diff-(tmax-tmin)>1e-5)
    {
      old_diff=tmax-tmin;
      
      float t=(tmax+tmin)/2;
      float z=T2z::getInstance().convert(t,dof);
      if (z>absz)
	tmax=t;
      else
	tmin=t;

      //cout << tmin << " " << tmax << endl;
    }

  if (z<0)
    tmin=-tmin;

  return tmin;
}

// }}}
// {{{ mygammapdf 

ReturnMatrix mygammapdf(const int npts, const float mult, const float delay, const float sigma)
{
  ColumnVector grot(npts);
  
  for(int i=0; i<npts; i++)
    grot(i+1)=i/mult;

  grot = gammapdf(grot.t(),delay,sigma*sigma).t();

  grot.Release();
  return grot;
}

// }}}
// {{{ estimate_X_heights 

ReturnMatrix estimate_X_heights(const Matrix real_X)
{
  ColumnVector real_X_heights(real_X.Ncols());

  for(int real_ev=0; real_ev<real_X.Ncols(); real_ev++)
    {
      // make sure that 0 is included in the X2 min:max range
      // this is so that (eg) all-1s EVs have height 1 (etc)
      float Xmax = Max( real_X.Column(real_ev+1).Maximum() , 0.0f );
      float Xmin = Min( real_X.Column(real_ev+1).Minimum() , 0.0f );

      real_X_heights(real_ev+1) = Xmax - Xmin;
    }

  real_X_heights.Release();
  return real_X_heights;
}

// }}}
// {{{ feat_svd 

ReturnMatrix feat_svd(const Matrix real_X)
{
  DiagonalMatrix eigenvals(1);
      
  if (real_X.Ncols()>1)
    {
      SVD(real_X,eigenvals);

      SortDescending(eigenvals);

      float inv_condition = eigenvals.Minimum() / eigenvals.Maximum();

      //cout << "eigenvals = " << eigenvals << endl;
      //cout << "inv_condition = " << inv_condition << endl;

      if (inv_condition<1e-7)
	{
	  cout << "Warning: at least one EV is (close to) a linear combination of the others. You probably need to alter your design.\n(Design matrix is rank deficient - ratio of min:max eigenvalues in SVD of matrix is " << inv_condition << ")\n" << endl;
	  //exit(0);
	}

    }

  eigenvals.Release();
  return eigenvals;
}

// }}}
// {{{ renorm kernel 

// void renorm_kernel(ColumnVector &X)
// {
//   double kernelsum=0;

//   for(int i=0; i<X.Nrows(); i++)
//     kernelsum += X(i+1); /* maybe this should be fabs(cX[i]), but looking at double-gamma output, I've not done this */

//   X /= kernelsum;
// }

// }}}
// {{{ orth_i_wrt_j 

// B -= B.A * A / A.A

void orth_i_wrt_j(Matrix &X, int i, int j)
{
  //cout << "orthogonalising real EV " << i << " wrt " << j << endl;

  float thedenom = ( X.Column(j).t() * X.Column(j) ).AsScalar();

  if (thedenom != 0)     // do nothing if Column(j) is constant
    X.Column(i) -= ( X.Column(i).t() * X.Column(j) ).AsScalar() * X.Column(j) / thedenom;
}

// }}}
// {{{ carry out the convolution 

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

// }}}
// {{{ find_line 

/* finds LAST matching entry in setup file */

char *find_line(char *filename, char *key, char *fl)
{
  FILE *fd;
  char *return_ptr=NULL, tmp_fl[10000];
  int  j;

  fd=fopen(filename,"rt");

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
      return (NULL);
    }
  else
    return return_ptr;
}

// }}}
// {{{ setup_font 

/* taken (I think) from pbmtext etc by Jef Poskanzer and George Phillips */

#define DEFAULTFONT_ROWS 155
#define DEFAULTFONT_COLS 112

#define FONT_WIDTH       7
#define FONT_HEIGHT      13

#define FONT_Y_SIZE     30

typedef struct { unsigned char font[DEFAULTFONT_ROWS][DEFAULTFONT_COLS];
                 int           char_row0[95], char_col0[95],
                               char_width, char_height; } FONT_DATA;

void error_exit(char *outkey )
{
  printf(outkey);
  return;
}

void setup_font(FONT_DATA *font_data)
{
    // {{{ Default Font 

/* The default font, packed in hex so this source file doesn't get huge.
   You can replace this with your own font using pbm_dumpfont().
*/

static unsigned long defaultfont_bits[DEFAULTFONT_ROWS][(DEFAULTFONT_COLS+31)/32] = {
    {0x00000000,0x20000c00,0x10000000,0x00000000},
    {0xc600a000,0x42000810,0x00000002,0x00000063},
    {0x6c00a000,0x45000810,0x00000002,0x00000036},
    {0x6c00a000,0x88800808,0xf2e1dee2,0x00000036},
    {0x54000000,0x80000800,0x11122442,0x0000002a},
    {0x54000001,0x00000800,0x11122442,0x0000002a},
    {0x54000001,0x00000800,0x11122282,0x0000002a},
    {0x44000102,0x00000800,0x11122382,0x00000022},
    {0xee000102,0x00000800,0x11e1e102,0x00000077},
    {0x00000204,0x00000800,0x11002102,0x00000000},
    {0x00000000,0x00000c00,0x11002102,0x00000000},
    {0x00000000,0x003f8000,0xe3807600,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x02000080,0x00040000,0x00120000,0x00000001},
    {0x04000082,0x828e1838,0x20210100,0x00000002},
    {0x04000082,0x82912448,0x20210100,0x00000002},
    {0x08000082,0x8fd01940,0x404087c2,0x00000004},
    {0x08000080,0x050c0622,0x00408102,0x00000004},
    {0x10000080,0x05061874,0x0040828f,0x00008008},
    {0x10000080,0x1f912688,0x00408002,0x00000008},
    {0x20000000,0x0a11098c,0x00408002,0x00000010},
    {0x20000080,0x0a0e0672,0x00210000,0x00000010},
    {0x40000000,0x00040000,0x00210000,0x00000020},
    {0x00000000,0x00000000,0x00120000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x004e0838,0x7023e1cf,0x00008000},
    {0x00000000,0x00913844,0x88620208,0x00008000},
    {0x08000000,0x00910844,0x08a20401,0x00000004},
    {0x10000000,0x01110844,0x08a20401,0x00000008},
    {0x20000000,0x01110808,0x3123c781,0x00000010},
    {0x400003e0,0x02110810,0x0a202441,0x00000020},
    {0x20000000,0x02110820,0x0bf02442,0x00000010},
    {0x10008000,0x04110844,0x88242442,0x00000008},
    {0x08008002,0x040e3e7c,0x7073c382,0x00000004},
    {0x00010000,0x08000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x0000e1c0,0x00000000,0x00000000,0x00000000},
    {0x00011220,0x00000000,0x70e38f87,0x00000000},
    {0x20011220,0x00020020,0x89108448,0x00008010},
    {0x10011220,0x00040010,0x09314448,0x00008008},
    {0x0800e221,0x02083e08,0x11514788,0x00000004},
    {0x040111e0,0x00100004,0x2153e448,0x00000002},
    {0x08011020,0x00083e08,0x213a2448,0x00008004},
    {0x10011040,0x02040010,0x01022448,0x00008008},
    {0x2000e381,0x02020020,0x20e77f87,0x00000010},
    {0x00000000,0x04000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x3803e7ef,0xc73bbe3d,0xdb863ce7,0x0000001c},
    {0x44011224,0x48910808,0x91036648,0x00008022},
    {0x4c011285,0x48910808,0xa1036648,0x00008026},
    {0x54011387,0x081f0808,0xc102a548,0x0000802a},
    {0x54011285,0x09910808,0xe102a548,0x0000802a},
    {0x4e011204,0x08910848,0x9112a4c8,0x00008027},
    {0x40011224,0x08910848,0x891224c8,0x00008020},
    {0x3803e7ef,0x073bbe31,0xcff77e47,0x0000001c},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000003,0x00000000},
    {0x0003e1cf,0x87bff7ef,0xdfbf77c2,0x00000000},
    {0x00013224,0x48a4a244,0x89122442,0x00000000},
    {0x00011224,0x4824a244,0xa8a14482,0x00000000},
    {0x00013227,0x8e04226c,0xa8414102,0x00000000},
    {0x0001e224,0x83842228,0xa8a08102,0x00000000},
    {0x00010224,0x40842228,0xd8a08242,0x00000000},
    {0x00010224,0x48843638,0x51108442,0x00000000},
    {0x0003c1ce,0x6f1f1c10,0x53b9c7c2,0x00000000},
    {0x00000060,0x00000000,0x00000002,0x00000000},
    {0x00000000,0x00000000,0x00000003,0x00000000},
    {0xfe000000,0x00000000,0x00000000,0x0000007f},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00010180,0x000000c0,0x003001c0,0x00000000},
    {0x08008081,0x00040040,0x00100200,0x00000004},
    {0x10008082,0x80040040,0x00100200,0x00000008},
    {0x10004084,0x40023c78,0x70f1c7c7,0x00004008},
    {0x10004080,0x00000244,0x89122208,0x00008008},
    {0x20002080,0x00001e44,0x8113e208,0x00008010},
    {0x10002080,0x00002244,0x81120208,0x00008008},
    {0x10001080,0x00002244,0x89122208,0x00008008},
    {0x10001080,0x00001db8,0x70e9c787,0x00008008},
    {0x10000880,0x00000000,0x00000000,0x00008008},
    {0x08000180,0x00000000,0x00000000,0x00008004},
    {0x00000000,0x1fc00000,0x00000007,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00030080,0x981c0000,0x00000000,0x00000000},
    {0x20010000,0x08040000,0x00000000,0x00000010},
    {0x10010000,0x08040000,0x00000000,0x00000008},
    {0x10016387,0x898474b8,0x72e1d5c7,0x00000008},
    {0x10019080,0x8a042a64,0x89122208,0x00008008},
    {0x08011080,0x8c042a44,0x89122207,0x00000004},
    {0x10011080,0x8a042a44,0x89122200,0x00008008},
    {0x10011080,0x89042a44,0x89122208,0x00008008},
    {0x1003bbe0,0x98dfebe6,0x71e1e787,0x00000008},
    {0x10000000,0x80000000,0x01002000,0x00000008},
    {0x20000000,0x80000000,0x01002000,0x00000010},
    {0x00000007,0x00000000,0x03807000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00008000,0x00000000,0x10410000,0x00000000},
    {0x00008000,0x00000000,0x20408000,0x00000000},
    {0x0001f66e,0xfdfbf77c,0x20408000,0x00000000},
    {0x24008224,0x488a2248,0x20408240,0x00000012},
    {0x54008224,0x4a842210,0x40404540,0x0000002a},
    {0x48008222,0x8a8a1420,0x20408480,0x00000024},
    {0x00008a23,0x85111c44,0x20408000,0x00000000},
    {0x000071d1,0x0531887c,0x20408000,0x00000000},
    {0x00000000,0x00000800,0x20408000,0x00000000},
    {0x00000000,0x00000800,0x10410000,0x00000000},
    {0x00000000,0x00003000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x00000000,0x00000000,0x00000000},
    {0x00000000,0x20000c00,0x10000000,0x00000000},
    {0xc600a000,0x42000810,0x00000002,0x00000063},
    {0x6c00a000,0x45000810,0x00000002,0x00000036},
    {0x6c00a000,0x88800808,0xf2e1dee2,0x00000036},
    {0x54000000,0x80000800,0x11122442,0x0000002a},
    {0x54000001,0x00000800,0x11122442,0x0000002a},
    {0x54000001,0x00000800,0x11122282,0x0000002a},
    {0x44000102,0x00000800,0x11122382,0x00000022},
    {0xee000102,0x00000800,0x11e1e102,0x00000077},
    {0x00000204,0x00000800,0x11002102,0x00000000},
    {0x00000000,0x00000c00,0x11002102,0x00000000},
    {0x00000000,0x003f8000,0xe3807600,0x00000000}
    };

// }}}
    unsigned char b;
    int rows, cols,
        scol, brow, d,
        bcol, ch;
    unsigned long l;

    // {{{ text explanation 

    /*
    ** This routine expects a font bitmap representing the following text:
    **
    ** (0,0)
    **    M ",/^_[`jpqy| M
    **
    **    /  !"#$%&'()*+ /
    **    < ,-./01234567 <
    **    > 89:;<=>?@ABC >
    **    @ DEFGHIJKLMNO @
    **    _ PQRSTUVWXYZ[ _
    **    { \]^_`abcdefg {
    **    } hijklmnopqrs }
    **    ~ tuvwxyz{|}~  ~
    **
    **    M ",/^_[`jpqy| M
    **
    ** The bitmap must be cropped exactly to the edges.
    **
    ** The dissection works by finding the first blank row and column; that
    ** gives the height and width of the maximum-sized character, which is
    ** not too useful.  But the distance from there to the opposite side is
    ** an integral multiple of the cell size, and that's what we need.  Then
    ** it's just a matter of filling in all the coordinates.
    */

// }}}
    // {{{ create font bits 

  for ( rows = 0; rows < DEFAULTFONT_ROWS; ++rows )
  {
    for ( cols = 0; cols < DEFAULTFONT_COLS; cols += 32 )
    {
      l = defaultfont_bits[rows][cols / 32];
      if (cols + 32 < DEFAULTFONT_COLS)
        ch = cols + 32;
      else
        ch = DEFAULTFONT_COLS;
      for ( scol = ch - 1; scol >= cols; --scol )
      {
        if ( l & 1 )
          font_data -> font[rows][scol] = 1;
        else
          font_data -> font[rows][scol] = 0;
        l >>= 1;
      }
    }
  }

// }}}
    // {{{ Find first blank row 

    for ( brow = 0; brow < DEFAULTFONT_ROWS / 6; ++brow )
	{
	b = font_data -> font[brow][0];
	for ( cols = 1; cols < DEFAULTFONT_COLS; ++ cols )
	    if ( font_data -> font[brow][cols] != b )
		goto nextrow;
	goto gotblankrow;
    nextrow: ;
	}

    error_exit("Couldn't find blank row in font\n");

// }}}
    // {{{ gotblankrow: 

gotblankrow:
    /* Find first blank col. */
    for ( bcol = 0; bcol < DEFAULTFONT_COLS / 8; ++bcol )
	{
	b = font_data -> font[0][bcol];
	for ( rows = 1; rows < DEFAULTFONT_ROWS; ++ rows )
	    if ( font_data -> font[rows][bcol] != b )
		goto nextcol;
	goto gotblankcol;
    nextcol: ;
	}
    error_exit("Couldn't find blank col in font.\n");

// }}}
    // {{{ gotblankcol: 

gotblankcol:
    /* Now compute character cell size. */
    d = DEFAULTFONT_ROWS - brow;
    font_data -> char_height = d / 11;
    if ( font_data -> char_height * 11 != d )
	error_exit("problem computing character cell height");
    d = DEFAULTFONT_COLS - bcol;
    font_data -> char_width = d / 15;
    if ( font_data -> char_width * 15 != d )
	error_exit("problem computing character cell width");

    /*printf("height=%d width=%d\n",font_data -> char_height,
           font_data -> char_width);*/

    /* Now fill in the 0,0 coords. */
    rows = font_data -> char_height * 2;
    cols = font_data -> char_width * 2;
    for ( ch = 0; ch < 95; ++ch )
	{
	font_data -> char_row0[ch] = rows;
	font_data -> char_col0[ch] = cols;
	cols += font_data -> char_width;
	if ( cols >= font_data -> char_width * 14 )
	    {
	    cols = font_data -> char_width * 2;
	    rows += font_data -> char_height;
	    }
	}

// }}}
}

// }}}
// {{{ write_string 

void write_string(unsigned char *in, int x, int y, char *the_string,
		  FONT_DATA *font_data, int colour, int x_size, int y_size)
{
  char string_bit;
  int  i, r, c, X, Y;

  for (i=0; the_string[i] != '\0'; i++)
  {
    string_bit = the_string[i] - ' ';
    for ( r = 0; r < font_data->char_height; ++r )
      for ( c = 0; c < font_data->char_width; ++c )
      {
        X = (i*font_data->char_width)+c+x;
        Y = r+y;
        if ( (X>-1) && (X<x_size) && (Y>-1) && (Y<y_size) &&
             (font_data->font[font_data->char_row0[(int)string_bit] + r][font_data->char_col0[(int)string_bit] + c]==1) )
          in[Y*x_size+X]=colour;
      }
  }
}

void write_string_rgb(unsigned char *r, unsigned char *g, unsigned char *b, 
		      int x, int y, char *the_string,
		      FONT_DATA *font_data, int cr, int cg, int cb,
		      int x_size, int y_size)
{
  write_string(r, x, y, the_string, font_data, cr, x_size, y_size);
  write_string(g, x, y, the_string, font_data, cg, x_size, y_size);
  write_string(b, x, y, the_string, font_data, cb, x_size, y_size);
}

// }}}

extern "C" __declspec(dllexport) int _stdcall feat_model(char *CmdLn)
{
  // {{{ variables 

   FILE   *ofp, *cofp;
   int    i, t, con, npts, ndelete, mnpts, level, *negpts, orig_evs, orig_ev, real_evs, real_ev, vox_evs, 
	   *orig_ev_nreal, *convolve_interaction, shape[10000], ncon, nftests, tempfiltyn, temphp_yn, templp_yn, f, *G=NULL;
   float  tr, mult, trmult, nltffwhm=0, tmpf, *triggers, maxconvwin=0;
   char   fn[10000], filename[10000], fl[10000], key[10000], the_string[10000], *FSLDIR;
   FONT_DATA *font_data = new FONT_DATA[1];
   DiagonalMatrix eigenvals;
   ColumnVector real_X_heights;

   int argc;
   char **argv;
  
   parser(CmdLn, argc, argv);

// }}}
  // {{{ read arguments and prepare variables 

if (argc<2)
{
  cout << "Usage: feat_model <design_name_root> [confound matrix text file]" << endl;
  freeparser(argc, argv);
  return(1);
}

Matrix motionparams(0,0);
if (argc==3)
  motionparams=remmean(read_ascii_matrix(argv[2]));
//cout << motionparams.Nrows() << "x" << motionparams.Ncols() << endl << motionparams << endl;

FSLDIR=getenv("FSLDIR");

sprintf(fn,"%s.fsf",argv[1]);

level     = atoi(find_line(fn, "fmri(level)", fl));
tr        = atof(find_line(fn, "fmri(tr)", fl));
npts      = atoi(find_line(fn, "fmri(npts)", fl));
ndelete   = atoi(find_line(fn, "fmri(ndelete)", fl));
orig_evs  = atoi(find_line(fn, "fmri(evs_orig)", fl))+(motionparams.Ncols()>0);
real_evs  = atoi(find_line(fn, "fmri(evs_real)", fl))+motionparams.Ncols();
vox_evs   = atoi(find_line(fn, "fmri(evs_vox)", fl));
ncon      = atoi(find_line(fn, "fmri(ncon_real)", fl));
nftests   = atoi(find_line(fn, "fmri(nftests_real)", fl));

npts    -= ndelete;     /* subtract off images which have already been deleted */
mult     = 1/DT;        /* points in high-res model per second */
trmult   = tr*mult;     /* points in high-res model per TR */
mnpts    = (int)(((float)(npts+10))*trmult+NEGSECS*mult); /* number of points in high-res model - add a few on for safety */

ColumnVector RE(ncon),                /* required effect (efficiency calculations) */
  real_CON_heights(ncon);

Matrix orig_X(mnpts,orig_evs);
orig_X=0;

ColumnVector orig_level(orig_evs); // records the difference between min and first point in an orig_X column

ColumnVector basisorth(orig_evs); // whether or not to orth basis functions wrt each other
basisorth=0;

Matrix real_X(npts,real_evs);
real_X=0;

Matrix C(real_evs,ncon);
Matrix F(nftests,ncon);

orig_ev_nreal = (int *)malloc(orig_evs*sizeof(int));

setup_font(font_data);

triggers = (float *)malloc(orig_evs*2*10*npts * sizeof(float)); /* factor of 10 for safety */
for(i=0;i<orig_evs*2*10*npts;i++)
  triggers[i]=-1e7;

negpts=(int*)calloc(orig_evs,sizeof(int));
convolve_interaction=(int*)calloc(orig_evs,sizeof(int));

float critical_t = z2t( atof(find_line(fn, "fmri(critical_z)", fl)) , MAX(npts-real_evs,1) );
float noise      = atof(find_line(fn, "fmri(noise)", fl));
float noisear    = atof(find_line(fn, "fmri(noisear)", fl));

/* setup prewhitening stuff */
SymmetricMatrix pwV(npts);
for(int i=1; i<=npts; i++)
  for(int j=1; j<=i; j++)
    pwV(i,j)=pow(noisear,(double)abs(i-j));
Matrix pwA = (Cholesky(pwV)).i(); // note newmat definition opposite of matlab
//cout << pwA * pwV * pwA.t() << endl; // should be identity

// }}}
  // {{{ read contrasts 

C=0;

sprintf(filename,"%s.con",argv[1]);
if ((cofp=fopen(filename,"wt"))==NULL)
{
  cout << "Can't open " << filename << " for writing" << endl;
  freeparser(argc, argv);
  return(1);
}

// modified from find_line
FILE *fd;
char tmp_fl[10000], *curconname;
fd=fopen(fn,"rt");
while ( fgets(tmp_fl, 1000, fd) != NULL ) {
  if (strncmp(tmp_fl,"set fmri(conname_real.",22)==0)
    {
      //tokenize tmp_vals...
      con=atoi(strtok(tmp_fl+22,") "));
      curconname=strtok(NULL,"\0\n\r");
      //still has newline, so don't need another one
      fprintf(cofp,"/ContrastName%d	%s",con,curconname);
    }
  if (strncmp(tmp_fl,"set fmri(con_real",17)==0)
    {
      //tokenize tmp_vals...
      con=atoi(strtok(tmp_fl+17,"."));// before period is contrastnum
      real_ev=atoi(strtok(NULL,") "));//between period and ") " is real_ev num
      C(real_ev,con)=atof(strtok(NULL," \0\n\t")); // rest is value
    }
  if (nftests>0)
    if (strncmp(tmp_fl,"set fmri(ftest_real",19)==0)
      {
	//tokenize tmp_vals...
	f=atoi(strtok(tmp_fl+19,"."));// before period is f-testnum
	con=atoi(strtok(NULL,") "));//between period and ") " is contrast num
	F(f,con)=atof(strtok(NULL," \0\n\t")); // rest is value
      }
 }
fclose(fd);

//check for 0 contrasts
for(con=1; con<=ncon; con++)
{
  int allzeros=1;

  for(real_ev=1; real_ev<=real_evs-motionparams.Ncols(); real_ev++)
    {
      if (C(real_ev,con)!=0) {
	allzeros=0;
	break;
      }
    }
  
  if (allzeros)
    {
      cout << "Contrast " << con << " is empty!" << endl;
      freeparser(argc, argv);
      return(1);
    }
}

// }}}
  // {{{ read and create design matrix 

if (level==1)
{
  // {{{ create basic shape 

for(orig_ev=0; orig_ev<orig_evs-(motionparams.Ncols()>0); orig_ev++)
  {
    orig_ev_nreal[orig_ev]=1;

    sprintf(key,"fmri(shape%d)",orig_ev+1); shape[orig_ev]=atoi(find_line(fn, key, fl));

    switch (shape[orig_ev])
      {
      case 0:
	// {{{ square wave 

{
  sprintf(key,"fmri(skip%d)",orig_ev+1);  float skip=atof(find_line(fn, key, fl))*mult;
  sprintf(key,"fmri(off%d)",orig_ev+1);   float off=atof(find_line(fn, key, fl))*mult;
  sprintf(key,"fmri(on%d)",orig_ev+1);    float on=atof(find_line(fn, key, fl))*mult;
  sprintf(key,"fmri(phase%d)",orig_ev+1); float phase=atof(find_line(fn, key, fl))*mult;
  sprintf(key,"fmri(stop%d)",orig_ev+1);  float stop=atof(find_line(fn, key, fl))*mult;

  if ( (stop<0) || (stop+skip>mnpts) ) stop=mnpts-skip;
  
  for(t=(int)skip;t<(int)(skip+stop);t++)
    { // do modulo maths in float not int - necessary for very short TR and block length
      if ( t+phase-skip - ((int)((t+phase-skip)/(off+on)))*(off+on) >= off )
	orig_X(t+1,orig_ev+1)=1.0;
      else
	orig_X(t+1,orig_ev+1)=0.0;
    }
}
break;

// {
//   sprintf(key,"fmri(skip%d)",orig_ev+1);  int skip=(int)(atof(find_line(fn, key, fl))*mult);
//   sprintf(key,"fmri(off%d)",orig_ev+1);   int off=(int)(atof(find_line(fn, key, fl))*mult);
//   sprintf(key,"fmri(on%d)",orig_ev+1);    int on=(int)(atof(find_line(fn, key, fl))*mult);
//   sprintf(key,"fmri(phase%d)",orig_ev+1); int phase=(int)(atof(find_line(fn, key, fl))*mult);
//   sprintf(key,"fmri(stop%d)",orig_ev+1);  int stop=(int)(atof(find_line(fn, key, fl))*mult);

//   if ( (stop<0) || (stop+skip>mnpts) ) stop=mnpts-skip;
  
//   for(t=skip;t<skip+stop;t++)
//     {
//       if ( (t+phase-skip)%(off+on) >= off )
// 	orig_X(t+1,orig_ev+1)=1.0;
//       else
// 	orig_X(t+1,orig_ev+1)=0.0;
//     }

// }
// break;

// }}}
      case 1:
	// {{{ sinusoid 

{
  sprintf(key,"fmri(skip%d)",orig_ev+1);   int skip=(int)(atof(find_line(fn, key, fl))*mult);
  sprintf(key,"fmri(period%d)",orig_ev+1); float period=atof(find_line(fn, key, fl)) * mult * 0.5 / M_PI;
  sprintf(key,"fmri(phase%d)",orig_ev+1);  float phase=atof(find_line(fn, key, fl)) * mult;
  sprintf(key,"fmri(stop%d)",orig_ev+1);   int stop=(int)(atof(find_line(fn, key, fl))*mult);

  if ( (stop<0) || (stop+skip>mnpts) ) stop=mnpts-skip;

  for(t=skip;t<skip+stop;t++)
    orig_X(t+1,orig_ev+1) = 0.5 * ( 1.0 + sin( M_PI + (t+phase-skip) / period ) );
}
break;

// }}}
      case 2:
	// {{{ custom single column 

{
  FILE *ifp2;

  sprintf(key,"fmri(custom%d)",orig_ev+1); strcpy(filename,(find_line(fn, key, fl)));

  if ((ifp2=fopen(filename,"rt"))==NULL)
    {
      cout << "Can't open " << filename << " for reading" << endl;
      freeparser(argc, argv);
      return(1);
    }

  for(i=0,t=0; t<npts; t++)
    {
      float tmpf;
	  
      if (fscanf(ifp2,"%f",&tmpf)!=1)
	{
	  cout << "Not enough data in " << filename << endl;
      freeparser(argc, argv);
	  return(1);
	}

      for(; i<trmult*(t+1); i++)
	orig_X(i+1,orig_ev+1)=tmpf;
    }

  fclose(ifp2);
}
break;

// }}}
      case 3:
	// {{{ custom 3 columns 

{
  FILE *ifp2;
  float start, stop, value;
  int   success=0;

  sprintf(key,"fmri(custom%d)",orig_ev+1); strcpy(filename,(find_line(fn, key, fl)));

  if ((ifp2=fopen(filename,"rt"))==NULL)
    {
      cout << "Can't open " << filename << " for reading" << endl;
      freeparser(argc, argv);
      return(1);
    }

  negpts[orig_ev] = (int)(NEGSECS*mult);

  while( fscanf(ifp2,"%f %f %f",&start,&stop,&value) == 3 )
    {
      start=start*mult+negpts[orig_ev]; stop=start+stop*mult;

      if ( (stop>0) && (stop-start<1) ) stop=start+1.1; /* if very short stim time make sure it happens */

      if (start<0)     start=0;    /* don't enter stuff before t=-negpts */
      if (start>mnpts) stop=0;     /* don't do anything */
      if (stop>mnpts)  stop=mnpts; /* don't overrun */

      for(t=(int)start; t<(int)stop; t++)
	{
	  orig_X(t+1,orig_ev+1)=value;
	  success=1;
	}
    }

  fclose(ifp2);

  if ( success==0 )
    {
      cout << "No valid\n[onset duration strength]\ntriplets found in " << filename << endl;
      freeparser(argc, argv);
      return(1);
    }

}
break;

// }}}
      case 4:
	// {{{ interactions

// complication; if the EVs used to feed into the interaction EV don't
// have the same convolution settings, then we cannot do this
// interaction before convolving the interaction - we must do the convolutions first....

{
  orig_X.Column(orig_ev+1)=1;
  int tmp, convEV=-1, other_convolves=-1;
  for(tmp=0; tmp<orig_ev; tmp++)
    {
      sprintf(key,"fmri(interactions%d.%d)",orig_ev+1,tmp+1);
      if ( atoi(find_line(fn, key, fl)) )
	{
	  float subtraction = orig_X(1,tmp+1) - orig_level(tmp+1);

	  sprintf(key,"fmri(convolve%d)",tmp+1); int convolve=atoi(find_line(fn, key, fl)); convEV=tmp;

	  if (other_convolves==-1)
	    other_convolves=convolve;
	  else
	    if ( other_convolves != convolve )
	      other_convolves=0;

	  ColumnVector tmpv(mnpts);
	  tmpv=0;
	  tmpv.Rows(1, mnpts-negpts[tmp]) = orig_X.SubMatrix(1+negpts[tmp], mnpts, tmp+1, tmp+1);

	  sprintf(key,"fmri(interactionsd%d.%d)",orig_ev+1,tmp+1);
	  int zero_how = atoi(find_line(fn, key, fl));
	  if ( zero_how == 1 )
	    subtraction = ( orig_X.Column(tmp+1).Maximum() + orig_X.Column(tmp+1).Minimum() ) / 2.0;
	  else if ( zero_how == 2 )
	    subtraction = orig_X.Column(tmp+1).Sum()/orig_X.Nrows();

	  // cout << "subtraction=" << subtraction <<endl;
	  orig_X.Column(orig_ev+1) = SP( orig_X.Column(orig_ev+1) , tmpv - subtraction );
	}
    }

  if (other_convolves>0)
    convolve_interaction[orig_ev]=convEV;
  //cout << orig_ev << " " << convolve_interaction[orig_ev] << endl;
}
break;
   case 9:
	// {{{ voxelwise input

{
  volume4D<float> ev_image;
  sprintf(key,"fmri(evs_vox_%d)",orig_ev+1);
  if ( read_volume4D(ev_image,find_line(fn, key, fl)) )
    cout << "Warning: voxelwise EV " << orig_ev+1 << " isn't readable" << endl;
  for(i=0,t=0; t<npts; t++)
    {
      float tmpf(ev_image[t].mean());
      for(; i<trmult*(t+1); i++)
	orig_X(i+1,orig_ev+1)=tmpf;
    }
 }
// }}}
      }

    //cout << orig_X << endl;

    orig_level(orig_ev+1) = orig_X(1,orig_ev+1) - orig_X.Column(orig_ev+1).Minimum();
    //cout << "EV" << orig_ev+1 << " min=" << orig_X.Column(orig_ev+1).Minimum() << " t[0]=" << orig_X(1,orig_ev+1) << " level=" << orig_level(orig_ev+1) << endl;

    if (shape[orig_ev]>1) // demean this orig EV
      orig_X.Column(orig_ev+1)=remmean(orig_X.Column(orig_ev+1));

    // {{{ triggers 

#define TRIGGER_THRESH 0.001

{
  float previous=0;
  int trigger_count=0;
  
  for(t=0;t<npts*trmult;t++)
    {
      if ( (orig_X(t+negpts[orig_ev]+1,orig_ev+1)>TRIGGER_THRESH) && (previous<=TRIGGER_THRESH) && (trigger_count<npts*2) )
	{
	  triggers[trigger_count*orig_evs+orig_ev]=t/trmult - 0.5; /* the 0.5 because sampling is halfway through TR */
	  trigger_count++;
	}

      if ( (orig_X(t+negpts[orig_ev]+1,orig_ev+1)<=TRIGGER_THRESH) && (previous>TRIGGER_THRESH) && (trigger_count<npts*2) )
	{
	  triggers[trigger_count*orig_evs+orig_ev]=t/trmult - 0.5; /* the 0.5 because sampling is halfway through TR */
	  trigger_count++;
	}

      previous=orig_X(t+negpts[orig_ev]+1,orig_ev+1);
    }

}

// }}}
  }
//Find first minimum timepoint for perfusion subtraction
 Matrix downsampledOriginalModel(npts,orig_evs);
 for(orig_ev=0; orig_ev<orig_evs-(motionparams.Ncols()>0); orig_ev++)
   do_resample(orig_X.Column(orig_ev+1),downsampledOriginalModel,orig_ev, trmult, negpts[orig_ev]);
 int minimumTimepoint=1;
 for (int row=1;row<=downsampledOriginalModel.Nrows();row++)
 {
   if ( downsampledOriginalModel.Row(row).Sum() < downsampledOriginalModel.Row(minimumTimepoint).Sum() )
     minimumTimepoint=row;
 }  
 ofstream outputFile((string(argv[1])+".min").c_str());
 if(outputFile.is_open()) {
   outputFile << minimumTimepoint;
   outputFile.close();
 }
// }}}
  // {{{ convolve and resample down in time 

for(orig_ev=real_ev=0; orig_ev<orig_evs-(motionparams.Ncols()>0); orig_ev++, real_ev++)
  if (shape[orig_ev]<10)
    {
      sprintf(key,"fmri(convolve%d)",orig_ev+1); int convolve=atoi(find_line(fn, key, fl));

      int c_orig_ev=orig_ev; // normally the same as orig_ev except for interactions
      if (shape[orig_ev] == 4) // i.e. interactions
	{
	  if (convolve_interaction[orig_ev]==0)
	    convolve=0;
	  else
	    {
	      c_orig_ev=convolve_interaction[orig_ev];
	      sprintf(key,"fmri(convolve%d)",c_orig_ev+1); convolve=atoi(find_line(fn, key, fl));
	    }
	}

      if (convolve>0)
	{
	  sprintf(key,"fmri(convolve_phase%d)",c_orig_ev+1); int convolve_phase=(int)(atof(find_line(fn, key, fl))*mult);

	  if ( (convolve>3) && (convolve<8) )
	    {
	      sprintf(key,"fmri(basisorth%d)",c_orig_ev+1);
	      basisorth(orig_ev+1)=atoi(find_line(fn, key, fl));
	    }

	  switch (convolve)
	    {
	    case 1:
	      // {{{ Gaussian 

{
  sprintf(key,"fmri(gausssigma%d)",c_orig_ev+1); float sigma=atof(find_line(fn, key, fl))*mult;
  sprintf(key,"fmri(gaussdelay%d)",c_orig_ev+1); float delay=atof(find_line(fn, key, fl))*mult;

  int fw = (int)(delay + sigma*5);
  maxconvwin=MAX(fw,maxconvwin);
  ColumnVector cX(fw);

  for(i=0; i<fw; i++)
    {
      float tmpf=(((float)i)-delay)/sigma;
      tmpf=exp(-0.5*tmpf*tmpf);
      if (i<1/DT) tmpf *= exp( -pow( 1-0.5*(((float)i)*DT) , 6.0 ) ); /* between 0 and 1s, window start of kernel */
      cX(i+1)=tmpf;
    }

  ColumnVector oX = do_convolve(orig_X.Column(orig_ev+1),cX,convolve_phase,1);
  do_resample(oX, real_X, real_ev, trmult, negpts[orig_ev]);
}
break;

// }}}
	    case 2:
	      // {{{ Gamma 

{
  sprintf(key,"fmri(gammasigma%d)",c_orig_ev+1); float sigma=atof(find_line(fn, key, fl));
  sprintf(key,"fmri(gammadelay%d)",c_orig_ev+1); float delay=atof(find_line(fn, key, fl));

  int fw = (int)((delay + sigma*5)*mult);
  maxconvwin=MAX(fw,maxconvwin);

  ColumnVector cX = mygammapdf(fw,mult,delay,sigma);

  ColumnVector oX=do_convolve(orig_X.Column(orig_ev+1),cX,convolve_phase,1);
  do_resample(oX, real_X, real_ev, trmult, negpts[orig_ev]);
}
break;

// }}}
	    case 3:
	      // {{{ double-gamma HRF 

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
break;

// }}}
	    case 4:
	      // {{{ Gamma basis functions 

#define WINDOW_FRAC 0.25

{
  sprintf(key,"fmri(basisfnum%d)",c_orig_ev+1);   int fnumber=atoi(find_line(fn, key, fl));
  sprintf(key,"fmri(basisfwidth%d)",c_orig_ev+1); int window=(int)(atof(find_line(fn, key, fl))*mult);

  orig_ev_nreal[orig_ev]+=fnumber-1;

  int fw = window;
  maxconvwin=MAX(fw,maxconvwin);

  float delay = fw/(2*mult);

  for(int fnum=1; fnum<=fnumber; fnum++)
    {
      float sigma = delay/2;

      ColumnVector cX = mygammapdf(fw,mult,delay,sigma);

      for(i=(int)((1-WINDOW_FRAC)*fw); i<fw; i++)
	cX(i+1) *= 1 - (i-(1-WINDOW_FRAC)*fw)/(WINDOW_FRAC*fw);

      ColumnVector oX=do_convolve(orig_X.Column(orig_ev+1),cX,convolve_phase,1);
      do_resample(oX, real_X, real_ev, trmult, negpts[orig_ev]);

      delay *= 0.5;
      real_ev++;
    }

  real_ev--;
}
break;

// }}}
	    case 5:
	      // {{{ Sine basis functions 

{
  sprintf(key,"fmri(basisfnum%d)",c_orig_ev+1);   int fnumber=atoi(find_line(fn, key, fl));
  sprintf(key,"fmri(basisfwidth%d)",c_orig_ev+1); int window=(int)(atof(find_line(fn, key, fl))*mult);

  orig_ev_nreal[orig_ev]+=fnumber-1;

  int fw = window;
  maxconvwin=MAX(fw,maxconvwin);
  ColumnVector cX(fw);

  for(int fnum=1; fnum<=fnumber; fnum++)
    {
      for(i=0; i<fw; i++)
	{
	  if (i<0) cout << "gcc is broken" << endl; // bizzarely, without this pointless line, the loop crashes....
	  cX(i+1)= sin(M_PI*i*fnum/fw);
	}

      ColumnVector oX=do_convolve(orig_X.Column(orig_ev+1),cX,convolve_phase,0);
      do_resample(oX, real_X, real_ev, trmult, negpts[orig_ev]);

      real_ev++;
    }

  real_ev--;
}
break;

// }}}
	    case 6:
	      // {{{ FIR basis functions 

{
  sprintf(key,"fmri(basisfnum%d)",c_orig_ev+1);   int fnumber=atoi(find_line(fn, key, fl));
  sprintf(key,"fmri(basisfwidth%d)",c_orig_ev+1); int window=(int)(atof(find_line(fn, key, fl))*mult);

  orig_ev_nreal[orig_ev]+=fnumber-1;

  int fw = window;
  maxconvwin=MAX(fw,maxconvwin);
  ColumnVector cX(fw);

  for(int fnum=0; fnum<fnumber; fnum++)
    {
      for(i=0; i<fw; i++)
	if ( (fw*fnum/fnumber<=i) && (i<fw*(fnum+1)/fnumber) )
	  cX(i+1)=1;
	else
	  cX(i+1)=0;

      ColumnVector oX=do_convolve(orig_X.Column(orig_ev+1),cX,convolve_phase,0);
      do_resample(oX, real_X, real_ev, trmult, negpts[orig_ev]);

      real_ev++;
    }

  real_ev--;
}
break;

// }}}
	    case 7:
	      // {{{ Custom basis functions 

{
  sprintf(key,"fmri(basisfnum%d)",c_orig_ev+1);   int fnumber=atoi(find_line(fn, key, fl));
  sprintf(key,"fmri(bfcustom%d)",c_orig_ev+1); char bfcustomname[10000]; strcpy(bfcustomname,find_line(fn, key, fl));

  Matrix icX(mnpts,fnumber);
  orig_ev_nreal[orig_ev]+=fnumber-1;
  
  FILE *ifp2;
  if ((ifp2=fopen(bfcustomname,"rt"))==NULL)
    {
      cout << "Can't open " << bfcustomname << " for reading" << endl;
      freeparser(argc, argv);
      return(1);
    }

  int fw, carryon=1;
  for(fw=0; carryon; fw++)
    for(int fnum=0; fnum<fnumber; fnum++)
      {
	float tmpval;
	if( fscanf(ifp2,"%f",&tmpval) != 1 )
	  carryon=0;
	icX(fw+1,fnum+1)=tmpval;
      }

  fw--;
  fclose(ifp2);
  maxconvwin=MAX(fw,maxconvwin);
  
  for(int fnum=0; fnum<fnumber; fnum++)
    {
      ColumnVector cX=icX.SubMatrix(1,fw,fnum+1,fnum+1);

      ColumnVector oX=do_convolve(orig_X.Column(orig_ev+1),cX,convolve_phase,0);
      do_resample(oX, real_X, real_ev, trmult, negpts[orig_ev]);

      real_ev++;
    }

  real_ev--;
}
break;

// }}}
	    }
	}
      else
	{
	  do_resample(orig_X.Column(orig_ev+1), real_X, real_ev, trmult, negpts[orig_ev]);
	  if ( shape[orig_ev]==1 )
	    // {{{ create sinusoid harmonics 

/* this is treated like a convolution to keep the structure simple -
thus the file key value searching is just a duplicate of that from the
principal sinusoid EV already created */

{
  int nharmonics;

  sprintf(key,"fmri(nharmonics%d)",c_orig_ev+1); nharmonics=atoi(find_line(fn, key, fl));

  if (nharmonics>0)
    {
      sprintf(key,"fmri(skip%d)",c_orig_ev+1);    int skip=(int)(atof(find_line(fn, key, fl))*mult);
      sprintf(key,"fmri(period%d)",c_orig_ev+1);  float period=atof(find_line(fn, key, fl)) * mult * 0.5 / M_PI;
      sprintf(key,"fmri(phase%d)",c_orig_ev+1);   float phase=atof(find_line(fn, key, fl)) * mult;
      sprintf(key,"fmri(stop%d)",c_orig_ev+1);    int stop=(int)(atof(find_line(fn, key, fl))*mult);

      if ( (stop<0) || (stop+skip>mnpts) ) stop=mnpts-skip;

      orig_ev_nreal[orig_ev]+=nharmonics;

      for(int harm=1; harm<=nharmonics; harm++)
	{
	  real_ev++;

	  for(t=skip;t<skip+stop;t++)
	    orig_X(t+1,orig_ev+1) = 0.5 * ( 1.0 + sin( M_PI + (harm+1)*(t+phase-skip) / period ) );
	
	  do_resample(orig_X.Column(orig_ev+1), real_X, real_ev, trmult, negpts[orig_ev]);
	}
    }
}

// }}}
	}
    }

// }}}
  // {{{ recomputed interactions

for(orig_ev=real_ev=0; orig_ev<orig_evs-(motionparams.Ncols()>0); real_ev+=orig_ev_nreal[orig_ev], orig_ev++)
  if ( (shape[orig_ev]==4) && (convolve_interaction[orig_ev]==0) )
    {
      real_X.Column(real_ev+1)=1;

      int tmp, tmp_real;
      for(tmp=tmp_real=0; tmp<orig_ev; tmp_real+=orig_ev_nreal[tmp], tmp++)
        {
          sprintf(key,"fmri(interactions%d.%d)",orig_ev+1,tmp+1);
          if ( atoi(find_line(fn, key, fl)) )
	    {
	      float subtraction = orig_X(1,tmp+1) - orig_level(tmp+1);

	      sprintf(key,"fmri(interactionsd%d.%d)",orig_ev+1,tmp+1);
	      int zero_how = atoi(find_line(fn, key, fl));
	      if ( zero_how == 1 )
		subtraction = ( real_X.Column(tmp_real+1).Maximum() + real_X.Column(tmp_real+1).Minimum() ) / 2.0;
	      else if ( zero_how == 2 )
		subtraction = real_X.Column(tmp_real+1).Sum()/real_X.Nrows();

	      // cout << "subtraction=" << subtraction <<endl;
	      real_X.Column(real_ev+1) = SP( real_X.Column(real_ev+1) , real_X.Column(tmp_real+1) - subtraction );
	    }
        }
    }

// }}}
  // {{{ add motion params to model if required 

if ( motionparams.Ncols() > 0 )
  {
    for(orig_ev=real_ev=0; orig_ev<orig_evs-1; real_ev+=orig_ev_nreal[orig_ev], orig_ev++);
    //cout << "inserting motion parameters starting at real EV " << real_ev+1 << endl;
    real_X.Columns(real_ev+1,real_ev+motionparams.Ncols()) = motionparams;
    orig_ev_nreal[orig_evs-1] = motionparams.Ncols();
    shape[orig_evs-1]=2;
  }

// }}}

  // two passes through the rest of the design setup:
  // first pass runs without temporal filtering, in order to get accurate peak-peak EV height estimation
  // second pass creates the final design, including temporal filtering
  Matrix tmp_real_X;
  int *tmp_orig_ev_nreal=(int*)malloc(sizeof(int)*orig_evs);
  for (int pass=0; pass<2; pass++)
    {
      if ( pass == 0 )
	{
	  tmp_real_X=real_X;
	  memcpy(tmp_orig_ev_nreal,orig_ev_nreal,sizeof(int)*orig_evs);
	} else {
	  real_X=tmp_real_X;
	  memcpy(orig_ev_nreal,tmp_orig_ev_nreal,sizeof(int)*orig_evs);
          // {{{ temporal filtering 

temphp_yn = atoi(find_line(fn, "fmri(temphp_yn)", fl));
templp_yn = atoi(find_line(fn, "fmri(templp_yn)", fl));
nltffwhm  = atof(find_line(fn, "fmri(paradigm_hp)", fl));

if ( (templp_yn) || (temphp_yn) )
  for(orig_ev=real_ev=0; orig_ev<orig_evs; orig_ev++)
  {
    if (orig_ev<orig_evs-(motionparams.Ncols()>0))
      {
	sprintf(key,"fmri(tempfilt_yn%d)",orig_ev+1);
	tempfiltyn=atoi(find_line(fn, key, fl));
      }
    else
      tempfiltyn=1;

    for(i=0; i<orig_ev_nreal[orig_ev]; i++, real_ev++)
      if (tempfiltyn)
	{
	  volume4D <float> im(1,1,1,npts);
	  double hp_sigma=-1, lp_sigma=-1;

	  if (templp_yn) lp_sigma=2.8/tr;
	  if (temphp_yn) hp_sigma=0.5*nltffwhm/tr;
	  /*printf("%f %f\n",hp_sigma,lp_sigma);*/

	  for(t=0;t<npts;t++)
	    im.value(0,0,0,t) = real_X(t+1,real_ev+1);

	  im=bandpass_temporal_filter(im, hp_sigma, lp_sigma);

	  for(t=0;t<npts;t++)
	    real_X(t+1,real_ev+1) = im.value(0,0,0,t);
	}
  }

// }}}
	}

      // {{{ demean 

real_X=remmean(real_X);

// }}}
      // {{{ orthogonalisation 

// basis function orthing (within orig_ev)
for(orig_ev=real_ev=0; orig_ev<orig_evs-(motionparams.Ncols()>0); real_ev+=orig_ev_nreal[orig_ev], orig_ev++)
  if (basisorth(orig_ev+1))
    {
      for(int tmp=1; tmp<orig_ev_nreal[orig_ev]; tmp++)
	{
	  Matrix tmp_orth = real_X.Columns(real_ev+1,real_ev+tmp);
	  real_X.Column(real_ev+tmp+1) = real_X.Column(real_ev+tmp+1) - tmp_orth*(pinv(tmp_orth)*real_X.Column(real_ev+tmp+1));
	}
    }


// old main orthogonalisation (WRONG)
// for(orig_ev=real_ev=0; orig_ev<orig_evs-(motionparams.Ncols()>0); real_ev+=orig_ev_nreal[orig_ev], orig_ev++)
// {
//   sprintf(key,"fmri(ortho%d.0)",orig_ev+1);
//   if (atoi(find_line(fn, key, fl)))
//     {
//       int tmp_orig_ev, tmp_real_ev;
//       for(tmp_orig_ev=tmp_real_ev=0; tmp_orig_ev<orig_evs-(motionparams.Ncols()>0); tmp_real_ev+=orig_ev_nreal[tmp_orig_ev], tmp_orig_ev++)
// 	if (tmp_orig_ev!=orig_ev)
// 	  {
// 	    sprintf(key,"fmri(ortho%d.%d)",orig_ev+1,tmp_orig_ev+1);
// 	    if (atoi(find_line(fn, key, fl)))
// 	      for(int tmp=0; tmp<MIN(orig_ev_nreal[orig_ev],orig_ev_nreal[tmp_orig_ev]); tmp++)
// 		orth_i_wrt_j(real_X,real_ev+tmp+1,tmp_real_ev+tmp+1);
// 	  }
//     }
// }



// main orthogonalisation
for(orig_ev=real_ev=0; orig_ev<orig_evs-(motionparams.Ncols()>0); real_ev+=orig_ev_nreal[orig_ev], orig_ev++)
{
  sprintf(key,"fmri(ortho%d.0)",orig_ev+1);
  if (atoi(find_line(fn, key, fl)))
    {
      Matrix tmp_orth(npts,0);
      int tmp_orig_ev, tmp_real_ev;
      for(tmp_orig_ev=tmp_real_ev=0; tmp_orig_ev<orig_evs-(motionparams.Ncols()>0); tmp_real_ev+=orig_ev_nreal[tmp_orig_ev], tmp_orig_ev++)
	if (tmp_orig_ev!=orig_ev)
	  {
	    sprintf(key,"fmri(ortho%d.%d)",orig_ev+1,tmp_orig_ev+1);
	    if (atoi(find_line(fn, key, fl)))  // i.e. should we orth orig_EV(orig_ev) wrt orig_EV(tmp_orig_ev)
	      for(int tmp=0; tmp<orig_ev_nreal[tmp_orig_ev]; tmp++)
		tmp_orth = tmp_orth | real_X.Column(tmp_real_ev+tmp+1);
	  }
      if (tmp_orth.Ncols()>0)
	for(int tmp=0; tmp<orig_ev_nreal[orig_ev]; tmp++)  // ORTH real_X.Column(real_ev+tmp+1) WRT tmp_orth		
	  real_X.Column(real_ev+tmp+1) = real_X.Column(real_ev+tmp+1) - tmp_orth*(pinv(tmp_orth)*real_X.Column(real_ev+tmp+1));
    }
}




// orthogonalisation wrt motion  (not used now - no point in any orthogonalisation!!)
// if (motionparams.Ncols()>0)
//   {
//     sprintf(key,"fmri(motionevs)");
//     int motionevs=atoi(find_line(fn, key, fl));

//     for(orig_ev=real_ev=0; orig_ev<orig_evs-1; real_ev+=orig_ev_nreal[orig_ev], orig_ev++);
//     //cout << real_ev << endl;

//     for (int i=0; i<real_ev; i++)
//       for (int j=real_ev; j<real_ev+motionparams.Ncols(); j++)
// 	{
// 	  if (motionevs==1)
// 	    orth_i_wrt_j(real_X,i+1,j+1);	  
// 	  else if (motionevs==2)
// 	    orth_i_wrt_j(real_X,j+1,i+1);	  
// 	}
//   }

// }}}
      // {{{ add temporal derivs 

for(orig_ev=real_ev=0; orig_ev<orig_evs-(motionparams.Ncols()>0); orig_ev++)
{
  sprintf(key,"fmri(deriv_yn%d)",orig_ev+1);

  if ( atoi(find_line(fn, key, fl)) )
    {
      orig_ev_nreal[orig_ev]++;

      // shift columns to the right to make space for tempderiv (if this isn't the right-most EV)
      if (orig_ev<orig_evs-1)
	{
	  //cout << "shifting real EVs " << real_ev+2 << ":" << real_evs-1 << " to " << real_ev+3 << ":" << real_evs << endl;
	  real_X.Columns(real_ev+3,real_evs) = real_X.Columns(real_ev+2,real_evs-1);
	}

      /* do end points explicitly */
      real_X(1,real_ev+2) = real_X(2,real_ev+1) - real_X(1,real_ev+1);
      real_X(npts,real_ev+2) = real_X(npts,real_ev+1) - real_X(npts-1,real_ev+1);

      for(t=1;t<npts-1;t++)
	real_X(t+1,real_ev+2) = (real_X(t+2,real_ev+1) - real_X(t,real_ev+1))*0.5;

      orth_i_wrt_j(real_X,real_ev+2,real_ev+1);	// well, Timmy wanted it but doesn't seem to have much effect.....
    }

  real_ev+=orig_ev_nreal[orig_ev];
}

// }}}
      // {{{ demean 

real_X=remmean(real_X);

// }}}

      if (pass==0)
	{
	  real_X_heights = estimate_X_heights(real_X);
	  // {{{ check rank of DM and get "real" contrast h2 heights (ie before HP filtering) 

// first do real rank deficiency test
// actually, no - we don't need this now we've switched to using pinv() below
//eigenvals = feat_svd(real_X);


// now do "meaningful" rank deficiency test

Matrix Q = pinv(real_X.t() * real_X);

for(con=1; con<=ncon; con++)
{
  ColumnVector contrast=C.Column(con);

  Matrix X2 = real_X * Q * contrast * pinv(contrast.t() * Q * contrast);

  float h2 = X2.Maximum() - X2.Minimum();
  if(fabs(h2)>1e6) h2=0; // try to catch dodgy heights, e.g. from empty EVs
  real_CON_heights(con)=h2;
}

// }}}
        }
      else
	{
          // {{{ output submodel group information 

{
  sprintf(filename,"%s.frf",argv[1]);
  if ((ofp=fopen(filename,"wt"))==NULL)
    {
      cout << "Can't open " << filename << " for writing" << endl;
      freeparser(argc, argv);
      return(1);
    }

  for(orig_ev=0; orig_ev<orig_evs; orig_ev++)
    for(int tmp=0; tmp<orig_ev_nreal[orig_ev]; tmp++)
      fprintf(ofp,"%d\n",orig_ev+1);

  fclose(ofp);
}

// }}}
	  // {{{ check rank of DM and do final contrast estimability test 

// first do real rank deficiency test
eigenvals = feat_svd(real_X);


// now do "meaningful" rank deficiency test

Matrix Q = pinv(real_X.t() * real_X);

for(con=1; con<=ncon; con++)
{
  float h2 = real_CON_heights(con); /* from first pass */

  ColumnVector dumbregressor = real_X * C.Column(con); /* just in order to test if it's empty */

  if ( dumbregressor.Maximum() - dumbregressor.Minimum() > 1e-10) /* i.e. not an 'empty' contrast, e.g. from empty EVs */
    {
      ColumnVector contrast=C.Column(con);

      Matrix X2 = real_X * Q * contrast * pinv(contrast.t() * Q * contrast);

      // whiten X2
      X2 = pwA*X2;

      //view whitening - don't leave this in!!
      //       real_X.Column(1)=X2;
      //       ncon=0;
      //       real_evs=orig_evs=1;
      //       orig_ev_nreal[0]=1;
      
      float D = h2 / sqrt( (X2.t() * X2).AsScalar() );
      RE(con) = critical_t * D * noise;
    }
  else
    {
      real_CON_heights(con)=0;
      C.Column(con)=0;
      RE(con) = 0;
    }
}

// }}}
        }
    }

} else {
  // {{{ group level design 

for(real_ev=0; real_ev<real_evs; real_ev++)
  {
    if (real_ev<real_evs-vox_evs)
      for(t=0;t<npts;t++)
	{
	  sprintf(key,"fmri(evg%d.%d)",t+1,real_ev+1);
	  real_X(t+1,real_ev+1)=atof(find_line(fn, key, fl));
	}
    else
      {
	volume4D<float> ev_image;
	sprintf(key,"fmri(evs_vox_%d)",1+real_ev-(real_evs-vox_evs));
	if ( read_volume4D(ev_image,find_line(fn, key, fl)) )
	  cout << "Warning: voxelwise EV " << 1+real_ev-vox_evs << " isn't readable" << endl;
	for(t=0;t<npts;t++)
	  real_X(t+1,real_ev+1)=ev_image[t].mean();
	  
      }
    orig_ev_nreal[real_ev]=1;
  }

// }}}
  // {{{ orthogonalisation 

for(real_ev=0; real_ev<real_evs-vox_evs; real_ev++)
  for(int tmp_real_ev=0; tmp_real_ev<real_evs-vox_evs; tmp_real_ev++)
    if (tmp_real_ev!=real_ev)
      {
	sprintf(key,"fmri(ortho%d.%d)",real_ev+1,tmp_real_ev+1);
	if (atoi(find_line(fn, key, fl)))
	  orth_i_wrt_j(real_X,real_ev+1,tmp_real_ev+1);
      }

// }}}
  real_X_heights = estimate_X_heights(real_X);
  // {{{ check rank of DM and do efficiency test 

// first do real rank deficiency test
eigenvals = feat_svd(real_X);
//eigenvals = feat_svd(real_X.Columns(1,real_evs-vox_evs));

// now do "meaningful" rank deficiency test

Matrix Q = pinv(real_X.t() * real_X);

for(con=1; con<=ncon; con++)
{
  ColumnVector contrast=C.Column(con);

  Matrix X2 = real_X * Q * contrast * pinv(contrast.t() * Q * contrast);

  // make sure that 0 is included in the X2 min:max range
  // this is so that (eg) all-1s EVs have height 1 (etc)
  float X2max = Max( X2.Maximum() , 0.0f );
  float X2min = Min( X2.Minimum() , 0.0f );

  float h2 = X2max - X2min;
  real_CON_heights(con)=h2;

  float D = h2 / sqrt( (X2.t() * X2).AsScalar() );

  RE(con) = critical_t * D * noise;
}

// }}}
}

// }}}
  // {{{ write matrix 

sprintf(filename,"%s.mat",argv[1]);
if ((ofp=fopen(filename,"wt"))==NULL)
{
  cout << "Can't open " << filename << " for writing" << endl;
  freeparser(argc, argv);
  return(1);
}

fprintf(ofp,"/NumWaves	%d\n",real_evs);
fprintf(ofp,"/NumPoints	%d\n",npts);

fprintf(ofp,"/PPheights	");
for(real_ev=1; real_ev<=real_evs; real_ev++)
  fprintf(ofp,"	%e",real_X_heights(real_ev));
fprintf(ofp,"\n");

fprintf(ofp,"\n/Matrix\n");

for(t=1; t<=npts; t++)
{
  for(real_ev=1; real_ev<=real_evs; real_ev++)
    fprintf(ofp,"%e	",real_X(t,real_ev));
  fprintf(ofp,"\n");
}

fclose(ofp);

// }}}
  // {{{ write triggers 

if (level==1)
{
  sprintf(filename,"%s.trg",argv[1]);
  if ((ofp=fopen(filename,"wt"))==NULL)
    {
      cout << "Can't open " << filename << " for writing" << endl;
      freeparser(argc, argv);
      return(1);
    }

  for(orig_ev=0; orig_ev<orig_evs; orig_ev++)
    {
      for(i=0;i<orig_ev_nreal[orig_ev];i++)
	{
	  if ( (shape[orig_ev]!=4) && (triggers[orig_evs+orig_ev]>-1e6) ) /* check for at least one trigger pair */
	    {
	      ColumnVector diffs_list(2*npts);
	      int j;

	      for(j=0;triggers[j*orig_evs+orig_ev]>-1e6;j+=2)
		{
		  fprintf(ofp,"%e ",triggers[j*orig_evs+orig_ev]);
		  if (triggers[(j+1)*orig_evs+orig_ev]>-1e6)
		    diffs_list(j/2+1)=triggers[(j+1)*orig_evs+orig_ev]-triggers[j*orig_evs+orig_ev];
		}

	      diffs_list=diffs_list.Rows(1,j/2);
	      fprintf(ofp,"%f\n",median(diffs_list) + maxconvwin/trmult);
	    }
	  else
	    fprintf(ofp,"0\n");
	}
    }

  fclose(ofp);
}

// }}}
  // {{{ write contrasts 

fprintf(cofp,"/NumWaves	%d\n",real_evs);
fprintf(cofp,"/NumContrasts	%d\n",ncon);

/* for contrasts, estimate "regressor height" as the mean of the heights of the relevant EVs */
fprintf(cofp,"/PPheights	");
for(con=1; con<=ncon; con++)
  fprintf(cofp,"	%e",real_CON_heights(con));
    
fprintf(cofp,"\n");

fprintf(cofp,"/RequiredEffect	");
for(con=1; con<=ncon; con++)
  fprintf(cofp,"	%.3f",RE(con));
fprintf(cofp,"\n");


fprintf(cofp,"\n/Matrix\n");
for(con=1; con<=ncon; con++)
{
  for(real_ev=1; real_ev<=real_evs; real_ev++)
    fprintf(cofp,"%e ",C(real_ev,con));
  fprintf(cofp,"\n");
}

fclose(cofp);

// }}}
  // {{{ ftests 

if (nftests>0)
{
  // new way -- add f-tests above with contrasts,
  // since that's loaded, come back and test here
  for(f=1; f<=nftests; f++)
    {
      Matrix Fmat(ncon,real_evs);
      int Fmat_rows=0;

      for(con=1; con<=ncon; con++)
	{
	  //	  sprintf(key,"fmri(ftest_real%d.%d)",f,con);
	  //F(f,con)=atoi(find_line(fn, key, fl));
	  if (F(f,con))
	    {
	      Fmat_rows++;
	      Fmat.Row(Fmat_rows)=(C.Column(con)).t();
	    }
	}

      Fmat=Fmat.SubMatrix(1,Fmat_rows,1,real_evs);

      // test that F(X'X)^-1F' is invertible, i.e. of full rank
	  if ((Fmat.Nrows() == 0) || (MISCMATHS::rank(Fmat*pinv(real_X.t()*real_X)*Fmat.t()) < Fmat.Nrows()))
	{
	  cout << "F-test " << f << " isn't valid - each included contrast cannot be a linear combination of the others." << endl;
      freeparser(argc, argv);
	  return(1);
	}
    }

  sprintf(filename,"%s.fts",argv[1]);
  if ((ofp=fopen(filename,"wt"))==NULL)
    {
      cout << "Can't open " << filename << " for writing" << endl;
      freeparser(argc, argv);
      return(1);
    }

  fprintf(ofp,"/NumWaves	%d\n",ncon);
  fprintf(ofp,"/NumContrasts	%d\n",nftests);
  fprintf(ofp,"\n/Matrix\n");

  for(f=1; f<=nftests; f++)
    {
      for(con=1; con<=ncon; con++)
	fprintf(ofp,"%d ",(int)F(f,con));
      fprintf(ofp,"\n");
    }

  fclose(ofp);
}

// }}}
  // {{{ second-level group memberships 

if (level==2)
{
  int maxG=0, isok, isnotzero;
  G=(int *)malloc(sizeof(int)*npts);

  for(t=0;t<npts;t++)
    {
      sprintf(key,"fmri(groupmem.%d)",t+1);
      G[t]=atoi(find_line(fn, key, fl));
      maxG=MAX(G[t],maxG);
    }

  sprintf(filename,"%s.grp",argv[1]);
  if ((ofp=fopen(filename,"wt"))==NULL)
    {
      cout << "Can't open " << filename << " for writing" << endl;
      freeparser(argc, argv);
      return(1);
    }

  /* if different group memberships check orthogonality of submatrices */
  if (maxG>1)
    {
      float *sub_X=(float *)malloc(sizeof(float)*npts*maxG);
      int *n_sub_X=(int *)malloc(sizeof(int)*maxG);
      for(real_ev=0; real_ev<real_evs; real_ev++)
	{
	  for(i=0;i<maxG;i++) n_sub_X[i]=0;
	  for(t=0;t<npts;t++)
	    sub_X[(G[t]-1)*npts+n_sub_X[G[t]-1]++] = real_X(t+1,real_ev+1);

	  isok=2;
	  for(i=0;i<maxG;i++)
	    {
	      isnotzero=0;
	      for(t=0;t<n_sub_X[i];t++)
		if (sub_X[i*npts+t]!=0)
		  isnotzero=1;
	      isok-=isnotzero;
	    }
	  if (isok<1)
	    printf("Warning - design matrix uses different groups (for different variances), but these do not contain \"separable\" EVs for the different groups (it is necessary that, for each EV, only one of the groups has non-zero values).\n");
	}
    }
  
  fprintf(ofp,"/NumWaves	1\n");
  fprintf(ofp,"/NumPoints	%d\n",npts);
  fprintf(ofp,"\n/Matrix\n");

  for(t=0;t<npts;t++)
    fprintf(ofp,"%d\n",G[t]);

  fclose(ofp);
}

// }}}
  // {{{ write covariance image 

{
  // {{{ setup vars 

FILE *ofp2;
unsigned char *r,*g,*b;
int border=5, mag=0, size=0, xsize=200, ysize, x, y, evx, evy;
float *cov=NULL, *covnorm=NULL;

if (real_evs>1)
  {
    cov = (float *)malloc(real_evs*real_evs*sizeof(float));
    covnorm = (float *)malloc(real_evs*sizeof(float));
    mag = MIN( MAX( 1 , 300/real_evs  ) , 50);
    size = real_evs*mag;
    xsize = size*2 + border*3;
  }

ysize = size   + border*2 + (level==1)*(border*(2+ncon) + FONT_HEIGHT*(ncon+1));

// }}}

  // {{{ malloc images and fill in background 

r=(unsigned char *)malloc(xsize*ysize);
g=(unsigned char *)malloc(xsize*ysize);
b=(unsigned char *)malloc(xsize*ysize);

memset((void *)r,(unsigned char)180,xsize*ysize);
memset((void *)g,(unsigned char)215,xsize*ysize);
memset((void *)b,(unsigned char)255,xsize*ysize);

// }}}

  if (real_evs>1)
  {						  
    // {{{ cov matrix orig 

for(evx=0; evx<real_evs; evx++)
  {
    for(evy=0; evy<real_evs; evy++)
      {
	cov[evy*real_evs+evx]= (real_X.Column(evx+1).t() * real_X.Column(evy+1)).AsScalar();
	//cout << cov[evy*real_evs+evx] << " ";
	if (evx==evy)
	  covnorm[evx]=cov[evx*real_evs+evx];
      }
    //cout << endl;
  }

for(evx=0; evx<real_evs; evx++)
  for(evy=0; evy<real_evs; evy++)
  {
    tmpf=covnorm[evy]*covnorm[evx];
    if (tmpf<=0)
      cov[evy*real_evs+evx]=0;
    else
      cov[evy*real_evs+evx]=abs(cov[evy*real_evs+evx])/sqrt(tmpf);
  }

for(evx=0; evx<real_evs; evx++)
  for(evy=0; evy<real_evs; evy++)
    for(y=0; y<mag; y++)
      for(x=0; x<mag; x++)
        r[(border+evy*mag+y)*xsize+border+evx*mag+x] =
        g[(border+evy*mag+y)*xsize+border+evx*mag+x] =
	  b[(border+evy*mag+y)*xsize+border+evx*mag+x] = (int)(255*cov[evy*real_evs+evx]);

// }}}
    // {{{ cov matrix SVD 

for(evx=0; evx<real_evs; evx++) for(evy=0; evy<real_evs; evy++) {
  int tmp=0;

  if (evx==evy)
    tmp=(int)(255*(eigenvals(evx+1)/eigenvals.Maximum()));

  for(y=0; y<mag; y++)
    for(x=0; x<mag; x++)
      r[(border+evy*mag+y)*xsize+2*border+size+evx*mag+x] =
	g[(border+evy*mag+y)*xsize+2*border+size+evx*mag+x] =
	b[(border+evy*mag+y)*xsize+2*border+size+evx*mag+x] = tmp;
}

for(evx=0; evx<real_evs; evx++) for(evy=0; evy<real_evs; evy++)
  if (evx==evy)
    {
      sprintf(the_string,"%.3f",eigenvals(evx+1)/eigenvals.Maximum());
      write_string_rgb(r, g, b, 
		       3*border+size+evx*mag + evx*(mag-2*border-5*FONT_WIDTH)/(real_evs-1),
		       border+evy*mag+mag/2-FONT_HEIGHT/2,
		       the_string, font_data, 255, 0, 0, xsize, ysize);
    }

// }}}
  }

  if (level==1)
    // {{{ contrast efficiencies 

{
  int yy;

  sprintf(the_string,"   Effect required (%%)");
  write_string_rgb(r, g, b, 
		   border,
		   size+3*border,
		   the_string, font_data, 0, 0, 0, xsize, ysize);

  for(yy=0; yy<ncon; yy++)
    {
      sprintf(the_string,"C%d %.3f",yy+1,RE(yy+1));
      
      int red_colour=0;
      /*if (RE(yy+1)>5)
	red_colour=255;*/

      write_string_rgb(r, g, b, 
		       border,
		       size+(4+yy)*border+(yy+1)*FONT_HEIGHT,
		       the_string, font_data, red_colour, 0, 0, xsize, ysize);
    }
}

// }}}

  // {{{ output image 

sprintf(filename,"%s_cov.ppm",argv[1]);
if ((ofp2=fopen(filename,"wb"))==NULL)
{
  cout << "Can't open " << filename << " for writing" << endl;
  freeparser(argc, argv);
  return(1);
}

fprintf(ofp2,"P6\n");
fprintf(ofp2,"%d %d\n",xsize,ysize);
fprintf(ofp2,"255\n");

for(y=0; y<ysize; y++)
     for(x=0; x<xsize; x++)
{
  fwrite(&r[y*xsize+x],1,1,ofp2);
  fwrite(&g[y*xsize+x],1,1,ofp2);
  fwrite(&b[y*xsize+x],1,1,ofp2);
}

fclose(ofp2);

   sprintf(filename,"wpng -q -overwrite %s_cov.ppm",argv[1]);
//   wpng(filename);
//system(filename);

// }}}
}

// }}}
  // {{{ write image preview (corrupts X[]) 

// {{{ setup vars 

FILE *ofp2;
unsigned char *r,*g,*b;
int xmag, ymag, border=5, temp=15, fsize=11, xsize, ysize, xx, yy, x, y, name_length=temp;

/* how much space for contrast names? */
for(yy=0; yy<ncon; yy++)
{
  sprintf(key,"fmri(conname_real.%d)",yy+1);
  strcpy(the_string,find_line(fn, key, fl));
  xx=(strlen(the_string)+4) * FONT_WIDTH;
  if (xx>name_length) name_length=xx;
}

xmag  = MIN( MAX( 1 , 600/real_evs  ) , 50);
ymag  = MIN( MAX( 1 , 400/npts ) , 20);

xsize=real_evs*xmag + border*(real_evs+3+nftests+(nftests>0)) + name_length + nftests*fsize;
ysize=npts*ymag + (ncon+1)*FONT_HEIGHT + border*(4+ncon);

// }}}
// {{{ reset X[] range (but don't change offset) 

for(real_ev=0; real_ev<real_evs; real_ev++)
  if ( real_X.Column(real_ev+1).MaximumAbsoluteValue() > 0 )
    real_X.Column(real_ev+1) /= real_X.Column(real_ev+1).MaximumAbsoluteValue();

// }}}
// {{{ malloc images and fill in background 

r=(unsigned char *)malloc(xsize*ysize);
g=(unsigned char *)malloc(xsize*ysize);
b=(unsigned char *)malloc(xsize*ysize);

memset((void *)r,(unsigned char)180,xsize*ysize);
memset((void *)g,(unsigned char)215,xsize*ysize);
memset((void *)b,(unsigned char)255,xsize*ysize);

// }}}
// {{{ time 

if (level==1)
for(yy=0; yy<npts; yy++)
  for(y=0; y<ymag; y++)
    for(x=0; x<temp; x++)
    {
      int intensity=64, rintensity;

      if (yy%10==0) intensity=255;

      rintensity=intensity;
      if (abs(yy-npts/2)<=(0.5*nltffwhm/tr))
	{
	  if (x>1&&x<temp-2)
	    rintensity=intensity=0;
	  if (x>3&&x<temp-4)
	    {
	      rintensity=255;
	      intensity=0;
	    }
	}

      r[(border+yy*ymag+y)*xsize+border+x]=rintensity;
      g[(border+yy*ymag+y)*xsize+border+x]=intensity;
      b[(border+yy*ymag+y)*xsize+border+x]=intensity;
    }

// }}}
// {{{ DM: grey 

  for(yy=0; yy<npts; yy++)
    for(xx=0; xx<real_evs; xx++)
      for(y=0; y<ymag; y++)
	for(x=0; x<xmag; x++)
	  {
	    float XX = real_X(yy+1,xx+1);

	    r[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+name_length]=(int)((XX * 100.0)+128);
	    g[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+name_length]=(int)((XX * 100.0)+128);
	    b[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+name_length]=(int)((XX * 100.0)+128);
	  }

// }}}
// {{{ DM: red 

  for(yy=0; yy<npts; yy++)
    for(xx=0; xx<real_evs; xx++)
      for(y=0; y<ymag; y++)
	{
	  float XX = real_X(yy+1,xx+1);
	  
	  x = (int)((XX*0.8+1.0)*xmag/2);

	  r[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+name_length]=255;
	  g[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+name_length]=0;
	  b[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+name_length]=0;
	  r[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+1+name_length]=255;
	  g[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+1+name_length]=0;
	  b[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+1+name_length]=0;
	  r[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-1+name_length]=255;
	  g[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-1+name_length]=0;
	  b[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-1+name_length]=0;
	  r[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+2+name_length]=0;
	  g[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+2+name_length]=0;
	  b[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x+2+name_length]=0;
	  r[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-2+name_length]=0;
	  g[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-2+name_length]=0;
	  b[(border+yy*ymag+y)*xsize+border*(xx+3)+xx*xmag+x-2+name_length]=0;
	}

// }}}
// {{{ EV names 

for(orig_ev=real_ev=0; orig_ev<orig_evs; orig_ev++)
  {
    if ( orig_ev < orig_evs-(motionparams.Ncols()>0) )
      {
	sprintf(key,"fmri(evtitle%d)",orig_ev+1);
	strcpy(the_string,find_line(fn, key, fl));
      }
    else
      sprintf(the_string,"%s","conf");

    for(i=0; i<orig_ev_nreal[orig_ev]; i++)
      {
	write_string_rgb(r, g, b, 
			 border*(real_ev+3)+real_ev*xmag+xmag/2+name_length-FONT_WIDTH*strlen(the_string)/2,
			 3*border+npts*ymag,
			 the_string, font_data, 150, 50, 50, xsize, ysize);
	real_ev++;
      }

  }

// }}}
// {{{ contrasts 

for(yy=0; yy<ncon; yy++)
{
  sprintf(the_string,"C%d",yy+1);
  write_string_rgb(r, g, b, 
		   border,
		   (4+yy)*border+(yy+1)*FONT_HEIGHT+npts*ymag,
		   the_string, font_data, 0, 0, 0, xsize, ysize);

  sprintf(key,"fmri(conname_real.%d)",yy+1);
  strcpy(the_string,find_line(fn, key, fl));
  write_string_rgb(r, g, b, 
		   border+4*FONT_WIDTH,
		   (4+yy)*border+(yy+1)*FONT_HEIGHT+npts*ymag,
		   the_string, font_data, 150, 50, 50, xsize, ysize);

  for(xx=0; xx<real_evs; xx++)
    {
      if ( C(xx+1,yy+1) == (int)C(xx+1,yy+1))
	sprintf(the_string,"%d",(int)C(xx+1,yy+1));
      else
	sprintf(the_string,"%.3f",C(xx+1,yy+1));
      
      write_string_rgb(r, g, b, 
		       border*(xx+3)+xx*xmag+xmag/2+name_length-FONT_WIDTH*strlen(the_string)/2,
		       (4+yy)*border+(yy+1)*FONT_HEIGHT+npts*ymag,
		       the_string, font_data, 0, 0, 0, xsize, ysize);
    }
}

// }}}
// {{{ ftests 

for(f=0; f<nftests; f++)
{
  write_string_rgb(r, g, b,
		   border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize),
		   3*border-FONT_HEIGHT+npts*ymag,
		   "F", font_data, 0, 0, 0, xsize, ysize);

  sprintf(the_string,"%d",f+1);
  write_string_rgb(r, g, b,
		   border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize),
		   3*border+npts*ymag,
		   the_string, font_data, 0, 0, 0, xsize, ysize);

  for(con=0; con<ncon; con++)
    {
      for(y=0; y<fsize; y++)
        for(x=0; x<fsize; x++)
	  r[ ((4+con)*border+(con+1)*FONT_HEIGHT+npts*ymag+y)*xsize +
	     border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize)+x ]=
            g[ ((4+con)*border+(con+1)*FONT_HEIGHT+npts*ymag+y)*xsize +
	       border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize)+x ]=
            b[ ((4+con)*border+(con+1)*FONT_HEIGHT+npts*ymag+y)*xsize +
	       border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize)+x ]=0;

      for(y=1; y<fsize-1; y++)
	for(x=1; x<fsize-1; x++)
	  if (F(f+1,con+1)==1)
	    r[ ((4+con)*border+(con+1)*FONT_HEIGHT+npts*ymag+y)*xsize +
	     border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize)+x ]=175;
	  else
	    {
	      r[ ((4+con)*border+(con+1)*FONT_HEIGHT+npts*ymag+y)*xsize +
		 border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize)+x ]=180;
	      g[ ((4+con)*border+(con+1)*FONT_HEIGHT+npts*ymag+y)*xsize +
		 border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize)+x ]=215;
	      b[ ((4+con)*border+(con+1)*FONT_HEIGHT+npts*ymag+y)*xsize +
		 border*(real_evs+4)+real_evs*xmag+name_length+f*(border+fsize)+x ]=255;
	    }
    }
}

// }}}
// {{{ second-level group memberships 

if (level==2)
  for(t=0;t<npts;t++)
    {
      sprintf(the_string,"%d",G[t]);
      write_string_rgb(r, g, b,
		       border,
		       border+t*ymag,
		       the_string, font_data, 0, 0, 0, xsize, ysize);
    }

// }}}
// {{{ output image 

sprintf(filename,"%s.ppm",argv[1]);
if ((ofp2=fopen(filename,"wb"))==NULL)
{
  cout << "Can't open " << filename << " for writing" << endl;
  freeparser(argc, argv);
  return(1);
}

fprintf(ofp2,"P6\n");
fprintf(ofp2,"%d %d\n",xsize,ysize);
fprintf(ofp2,"255\n");

for(y=0; y<ysize; y++)
     for(x=0; x<xsize; x++)
{
  fwrite(&r[y*xsize+x],1,1,ofp2);
  fwrite(&g[y*xsize+x],1,1,ofp2);
  fwrite(&b[y*xsize+x],1,1,ofp2);
}

fclose(ofp2);
sprintf(filename,"wpng -q -overwrite  %s.ppm ",argv[1]);
//wpng(filename);

//sprintf(filename,"%s/bin/wpng -q -overwrite  %s.ppm ",getenv("FSLDIR"),argv[1]);
//system(filename);

// }}}

// }}}

  freeparser(argc, argv);
  return(0);
}
}
