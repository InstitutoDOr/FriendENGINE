//     fslnvols.cc   REALLY advanced program for counting the volumes in a 4D file
//     Stephen Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 1999-2005 University of Oxford  
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

#include "newimage/newimageall.h"

using namespace NEWIMAGE;

namespace fslnvols {
#include "newimage/fmribmain.h"
void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: fslnvols <input>" << endl;
}


template <class T>
int fmrib_main(int argc, char *argv[])
{
  FSLIO *IP1;
  IP1 = NewFslOpen(string(argv[1]), "r");
  if (IP1==0) throw Exception(("Failed to read volume "+string(argv[1])).c_str()); 
  short sx=0,sy=0,sz=0,st=0;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  FslClose(IP1);  
  cout << st << endl;
  return 0;
}



int main(int argc,char *argv[])
{
  string progname=argv[0];
  if (argc != 2)
  {
    print_usage(progname);
    return 1;
  }
  if (!FslFileExists(argv[1])) 
  { 
    cout << "0" << endl;
    return 0; // the intended usage is to return "0" and not show an error
  }
  else 
  { 
    string iname=string(argv[1]);
    return call_fmrib_main(dtype(iname),argc,argv); 
  }
}

}