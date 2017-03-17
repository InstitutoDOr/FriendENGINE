/*  CCOPYRIGHT  */
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

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

namespace cutoffcalc {
string title="Cut-off calcualtor (Version 1.1)\nCopyright(c) 2009, University of Oxford (Mark Jenkinson and Matthew Webster)\nCalculates the minimal period for the highpass filter that still preserves a specified amount of variance in all the design matrix regressors.";
string examples="cutoffcalc [options] -i design.mat";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false,
		     string("switch on diagnostic messages"),
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<float> varthreshold(string("-t,--thresh"), 0.9,
		  string("Threshold for retained variance (default=0.9)"),
		  false, requires_argument);
Option<float> tr(string("--tr"), 3,
		  string("Time between successive data points (default=3.0s)"),
		  false, requires_argument);
Option<float> limit(string("--limit"), 90.0,
		  string("Lower limit on period due to autocorr estimation (default=90s)"),
		  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input design matrix"),
		  true, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////
double min_std_ratio(double sigma, volume4D<float>& vorig, const volume<float>& std0)
{
  static volume4D<float> vtmp;
  static volume<float> std1;
  double stdratio(1.0);
  vtmp=bandpass_temporal_filter(vorig,sigma,0.0);
  std1 = stddevvol(vtmp);
  float maxstd=std0.max();
  for (int z=0; z<=std0.maxz(); z++) { 
    for (int y=0; y<=std0.maxy(); y++) { 
      for (int x=0; x<=std0.maxx(); x++) { 
	if (std0(x,y,z)>1e-12*maxstd) {  // in case there are zero regressors
	  std1(x,y,z) = std1(x,y,z) / std0(x,y,z);
	} else {
	  std1(x,y,z) = 1.0;  // always passes test
	}
      }
    }
  }
  stdratio=std1.min();
  return stdratio; 
}

int do_work(int argc, char* argv[])
{
 Matrix dm(read_vest(inname.value()));
 remmean(dm);
 volume4D<float> vorig;
 volume<float> std0, mask(dm.Ncols(),1,1);
 mask=1.0;
 vorig.setmatrix(dm,mask);
 std0 = stddevvol(vorig);
 float sig2sec = tr.value()*2.0;
 float lowersig = limit.value()/sig2sec;  // convert secs to sigma
 float uppersig = 3*dm.Nrows();
 float sigma(uppersig), stdratio(1.0);
 float usig(uppersig), lsig(lowersig), hsig;
 // sanity check at ends
 if (min_std_ratio(usig,vorig,std0)<=varthreshold.value()) {
   if (verbose.value()) { cout << "Failed to meet threshold criterion: extreme period still removes too much variance" << endl; } 
   cout << usig*sig2sec << endl;
   return 0;
 }
 if (min_std_ratio(lsig,vorig,std0)>varthreshold.value()) {
   if (verbose.value()) { cout << "Using lowest period in range" << endl; }
   cout << lsig*sig2sec << endl;
   return 0;
 }
 do {
   hsig=MISCMATHS::round((lsig+usig)/2);  // force hsig to equal usig when lsig and usig differ by 1
   stdratio=min_std_ratio(hsig,vorig,std0);
   if (stdratio>varthreshold.value()) { usig=hsig; }
   else { lsig=hsig; }
 } while ((usig-lsig)>1);
 
 if (verbose.value()) { cout << "Sigma bounds are " << lsig << " and " << usig << endl; }
 if (verbose.value()) { cout << "stdratios are " << min_std_ratio(lsig,vorig,std0) << " and " << min_std_ratio(usig,vorig,std0) << endl; }

 sigma=usig;

 if (verbose.value()) { cout << "FEAT highpass filter value (in seconds) should be: " << endl; }
 cout << sigma*sig2sec << endl;
 volume4D<float> vtmp;
 if ( verbose.value() ) { 
   vtmp=bandpass_temporal_filter(vorig,sigma,0);
   float baseSum=(dm-vtmp.matrix(mask)).SumSquare();
   Matrix base=vtmp.matrix(mask);
   for (float factor=0.6;factor<=1.6;factor+=0.1)
   {
     vtmp=bandpass_temporal_filter(vorig,sigma*factor,0);
     if (verbose.value()) {
       //cout << factor << "\t" << ((dm-vtmp.matrix(mask)).SumSquare())/baseSum << endl;
       cout << factor << "\t" << ((base-vtmp.matrix(mask)).SumSquare())/baseSum << endl;
     }
   }
 }

 return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

 Tracer tr0("main");
 OptionParser options(title, examples);

 try {
   // must include all wanted options here (the order determines how
   //  the help message is printed)
   options.add(inname);
   options.add(varthreshold);
   options.add(tr);
   options.add(limit);
   options.add(verbose);
   options.add(help);

   nonoptarg = options.parse_command_line(argc, argv);

   // line below stops the program if the help was requested or
   //  a compulsory option was not set
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

 return do_work(argc,argv);
}
}