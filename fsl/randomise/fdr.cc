/*  fdr.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2004-2006 University of Oxford  */

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

// Calculates the p-threshold or minimum q rate image from an
//  input 3D probability image (values of 0 for p are ignored)

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include <vector>
#include <algorithm>
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

namespace fdr {
string title="fdr (Version 1.2)\nCopyright(c) 2004-2006, University of Oxford (Mark Jenkinson)";
string examples="fdr -i <pvalimage> [options]\ne.g.  fdr -i <pvalimage> -m <maskimage> -q 0.05\n      fdr -i <pvalimage> -o <qrateimage>\nNote: if a mask is not specified, voxels where p>.9999 are ignored.";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> conservativetest(string("--conservative"), false,
			      string("use conservative FDR correction factor"),
			      false, no_argument);
Option<float>  qthresh(string("-q"),0,
		  string("q-value (FDR) threshold"),
		  false, requires_argument);
Option<string> ordername(string("--order"), string(""),
		      string("~\toutput image of order values"),
		      false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		      string("input p-value filename"),
		      true, requires_argument);
Option<string> mask(string("-m"), string(""),
		      string("mask filename"),
		      false, requires_argument);
Option<string> qoutname(string("-o"), string(""),
		       string("output q-rate filename"),
		       false, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

vector<int> get_sortindex(const Matrix& vals)
{
  // return the mapping of old indices to new indices in the
  //   new *ascending* sort of vals
  int length=vals.Nrows();
  vector<pair<double, int> > sortlist(length);
  for (int n=0; n<length; n++) {
    sortlist[n] = pair<double, int>((double) vals(n+1,1),n+1);
  }
  sort(sortlist.begin(),sortlist.end());  // O(N.log(N))
  vector<int> idx(length);
  for (int n=0; n<length; n++) {
    idx[sortlist[n].second-1] = n+1;
  }
  return idx;
}


////////////////////////////////////////////////////////////////////////////

int save_as_image(const string& filename, const volume<float>& mask, 
		  const Matrix& valmat)
{
    // put values back into volume format
    if (verbose.value()) { cerr << "Saving results to " << filename << endl; }
    volume4D<float> outvals;
    outvals.addvolume(mask);
    outvals.setmatrix(valmat.t(),mask);
    return save_volume4D(outvals,filename);
}


int do_work(int argc, char* argv[], int nonoptarg) 
{
  volume4D<float> pimg; // reading in 4D in order to be able to use the newimage<->newmat functions - only actually use first timepoint
  read_volume4D(pimg,inname.value());
  if (verbose.value()) print_info(pimg,"p-value image");

  volume<float> vmask;
  if (mask.set())
    {
      read_volume(vmask,mask.value());
      vmask.binarise(0.0001);
    }
  else
    {
      vmask = pimg[0];
      vmask.binarise(0.9999);
      vmask = 1.0f-vmask;
    }

  Matrix pmat;
  pmat = pimg.matrix(vmask);
  pmat = pmat.t();
  int Ntot = pmat.Nrows();

  // calculate FDR threshold required to make each voxel significant
  // FDR formula is: p = n*q / (N * C)  
  //   where n=order index, N=total number of p values, 
  //         C=1 for the simple case, and 
  //         C=1/1 + 1/2 + 1/3 + ... + 1/N for the most general correlation
  // We use the inverse formula: q_{min} = N*C*p / n

  if (verbose.value()) { cerr << "Calculating FDR values" << endl; } 
  float C=1.0;
  if (conservativetest.value()) {
    for (int n=2; n<=Ntot; n++) { C+=1.0/((double) n); }
  }

  vector<int> norder = get_sortindex(pmat);

  // output the appropriate p-threshold, if requested
  if (qthresh.set()) {
    float qthr = qthresh.value();
    float qfac = qthr / ( C * (float) Ntot );
    float pthresh = 0.0;
    for (int j=1; j<=Ntot; j++) {
      if ( (pmat(j,1) > pthresh) && 
	   ( pmat(j,1) < qfac * (float) norder[j-1] ) )
      {
	 if (verbose.value()) { cout << "p = " << pmat(j,1) << " , n = " 
	         << norder[j-1] << " , qfac = " << qfac << endl; }
         pthresh = pmat(j,1);
      }
    }
    cout << "Probability Threshold is: " << endl << pthresh << endl;
  }
  

  // output the q (fdr) image, if requested
  if (qoutname.set()) {
    Matrix qmat = pmat;
    for (int j=1; j<=Ntot; j++) {
      qmat(j,1) = pmat(j,1) * Ntot * C / norder[j-1];
    }
    
    // save the FDR (q_min) results
    save_as_image(qoutname.value(),vmask,qmat);
  }

  
  // save the order values, if requested
  if (ordername.set()) {
    Matrix ordermat(Ntot,1);
    for (int j=1; j<=Ntot; j++) { ordermat(j,1) = norder[j-1]; }
    save_as_image(ordername.value(),vmask,ordermat);
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(inname);
    options.add(mask);
    options.add(qthresh);
    options.add(qoutname);
    options.add(ordername);
    options.add(conservativetest);
    options.add(verbose);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions

  return do_work(argc,argv,nonoptarg);
}

}