/*  HalfcosbasisOptions.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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

#if !defined(HalfcosbasisOptions_h)
#define HalfcosbasisOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace Halfcosbasis {

class HalfcosbasisOptions {
 public:
  static HalfcosbasisOptions& getInstance();
  ~HalfcosbasisOptions() { delete gopt; }
  
  Option<bool> verbose;
  Option<int> debuglevel;
  Option<bool> timingon;
  Option<bool> help;
  Option<string> logdir;
  Option<string> halfcosparamrangesfile;
  Option<int> nhrfsamps;
  Option<int> nbfs;
  Option<int> nsecs;
  Option<float> res;

  void parse_command_line(int argc, char** argv, Log& logger);
  
 private:
  HalfcosbasisOptions();  
  const HalfcosbasisOptions& operator=(HalfcosbasisOptions&);
  HalfcosbasisOptions(HalfcosbasisOptions&);

  OptionParser options; 
      
  static HalfcosbasisOptions* gopt;
  
};

 inline HalfcosbasisOptions& HalfcosbasisOptions::getInstance(){
   if(gopt == NULL)
     gopt = new HalfcosbasisOptions();
   
   return *gopt;
 }

 inline HalfcosbasisOptions::HalfcosbasisOptions() :
   verbose(string("-V,--verbose"), false, 
	   string("switch on diagnostic messages"), 
	   false, no_argument),
   debuglevel(string("--db,--debug,--debuglevel"), 0,
		       string("set debug level"), 
		       false, requires_argument),
   timingon(string("--to,--timingon"), false, 
		       string("turn timing on"), 
		       false, no_argument),
   help(string("-h,--help"), false,
		    string("display this message"),
		    false, no_argument),
   logdir(string("-l,--ld,--logdir"), string("logdir"),
	  string("log directory"),
	  false, requires_argument),  
   halfcosparamrangesfile(string("--hcprf,--hf"), string(""),
	      string("Half cosine HRF parameter ranges file"),
	      true, requires_argument),
   nhrfsamps(string("--nhs"), 1000,
	  string("Num of HRF samples to use (default is 1000)"),
	  false, requires_argument), 
   nbfs(string("--nbfs"), 3,
	  string("Num of HRF basis functions to use (default is 3)"),
	  false, requires_argument),
   nsecs(string("--ns,--nsecs"), 40,
	  string("Num of secs (default is 40)"),
	  false, requires_argument), 
   res(string("--res"), 0.05,
	  string("Temporal resolution (default is 0.05)"),
	  false, requires_argument), 
   options("halfcosbasis","halfcosbasis --help")
 {
     try {
       options.add(verbose);
       options.add(debuglevel);
       options.add(timingon);
       options.add(logdir);
       options.add(halfcosparamrangesfile);
       options.add(nhrfsamps);
       options.add(nbfs);
       options.add(nsecs);
       options.add(res);
       options.add(help);
     }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
     
   }
}

#endif



