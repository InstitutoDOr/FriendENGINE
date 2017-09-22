/*  fslfft.cc

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2003-2007 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
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
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    Innovation@innovation.ox.ac.uk quoting reference DE/9564. */


// Calculates and prints the nth cog of an image

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "parser.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;


namespace fslfft {
void print_usage(const string& progname) {
  cerr << "Usage: " << progname << " <input volume> <output volume> [-inv]" 
       << endl << endl;
  cerr << "e.g.   " << progname << " invol outvol" 
       << endl;
}



extern "C" __declspec(dllexport) int _stdcall fslfft(char *CmdLn)
{
  int argc;
  char **argv;
  
  parser(CmdLn, argc, argv);

  Tracer tr("main");

  string progname=argv[0];
  if (argc<3) { 
    print_usage(progname);
    return -1; 
  }

  bool inverse=false;
  if (argc>=4) {
    string optarg=argv[3];
    if (optarg=="-inv") {
      inverse=true;
    } else {
      cerr << "Unrecognised option: " << argv[3] << endl;
      return 1;
    }
  }

  string inname=string(argv[1]), outname=string(argv[2]);


  complexvolume vin, vc;
  read_complexvolume(vin,inname);
  vc = vin;
  if (inverse) {
    ifft2(vc);
  } else {
    fft2(vc);
  }
  save_complexvolume(vc,outname);

  freeparser(argc, argv);
  return 0;
}
}
