/*  FLAME - FMRIB's Local Analysis of Mixed Effects

    Mark Woolrich, Tim Behrens - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH
#include "newmatap.h"
#include "newmatio.h"
#include <string>
#include <math.h>
#include "utils/log.h"
#include "gsmanager.h"
#include "gsoptions.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscprob.h"
#include "stdlib.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace Gs;
using namespace MISCMATHS;

int main(int argc, char *argv[])
{
  try{  

//     for(int i=1; i < 10; i++)
//       {
// 	OUT(rand()/double(RAND_MAX));
//         OUT(unifrnd());
//         OUT(normrnd());
//       }
//     exit(0);


    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    
    // parse command line - will output arguments to logfile
    GsOptions& opts = GsOptions::getInstance();
    opts.parse_command_line(argc, argv, logger);

    srand(GsOptions::getInstance().seed.value());
   
    if(opts.debuglevel.value()==1)
      Tracer_Plus::setrunningstackon();

    if(opts.timingon.value())
      Tracer_Plus::settimingon();


//     ColumnVector storen(10000);
//     ColumnVector store(10000);
//     for(int i=1; i<=10000; i++)
//       {
// 	storen(i) = normrnd().AsScalar();
// 	store(i) = (unifrnd().AsScalar()+0.5);
//       }

//     write_ascii_matrix(storen,"storen");
//     write_ascii_matrix(store,"store");

    Gsmanager gsmanager;
    cout << "Setting up:" << endl;
    gsmanager.setup();
    cout << "Running:" << endl;
    gsmanager.run();
    cout << "Saving results" << endl;
    gsmanager.save();

    if(opts.timingon.value())
      Tracer_Plus::dump_times(logger.getDir());

    cout << endl << "Log directory was: " << logger.getDir() << endl;
  }
  catch(Exception& e) 
    {
      cerr << endl << e.what() << endl;
      return 1;
    }
  catch(X_OptionError& e) 
    {
      cerr << endl << e.what() << endl;
      return 1;
    }

  return 0;
}












