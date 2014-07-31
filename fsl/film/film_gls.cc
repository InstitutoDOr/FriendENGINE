/*  film_gls.cc

    Mark Woolrich and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2008 University of Oxford  */

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
#define WANT_STREAM
#define WANT_MATH

#include "newimage/newimageall.h"
#include "utils/log.h"
#include "AutoCorrEstimator.h"
#include "paradigm.h"
#include "FilmGlsOptions.h"
#include "glimGls.h"
#include "parser.h"
#include <vector>
#include <string>

using namespace NEWMAT;
using namespace FILM;
using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

extern "C" __declspec(dllexport) int _stdcall film_gls(char *CmdLn)
{
  int argc;
  char **argv;

  parser(CmdLn, argc, argv);

  try{
    
    // Setup logging:
    Log& logger = LogSingleton::getInstance();

    // parse command line
    FilmGlsOptions& globalopts = FilmGlsOptions::getInstance();
    globalopts.parse_command_line(argc, argv, logger);

    // load data
    volume4D<float> input_data;
    read_volume4D(input_data,globalopts.inputfname);
    int sizeTS(input_data.tsize());

    volume4D<float> reference;
    reference=input_data[int(sizeTS/2)-1];
    copybasicproperties(input_data,reference);

    volume<float> mask(meanvol(input_data));
    volume<float> variance(variancevol(input_data));
    input_data-=mask;

    mask.binarise(globalopts.thresh,mask.max()+1,exclusive);
    variance.binarise(1e-10,variance.max()+1,exclusive); //variance mask needed if thresh is -ve to remove background voxels (0 variance)
    mask*=variance; //convolved mask ensures that only super-threshold non-background voxels pass

    Matrix datam(input_data.matrix(mask));
   
    int numTS = datam.Ncols();
    ColumnVector epivol = reference.matrix(mask).t();
    
    // Load paradigm:
    Paradigm paradigm;
    if(!globalopts.ac_only)
      paradigm.load(globalopts.paradigmfname, "", "", false, sizeTS);
    else
      paradigm.setDesignMatrix(sizeTS); // set design matrix to be one ev with all ones:

    if(globalopts.verbose)
      write_vest(logger.appendDir("Gc"), paradigm.getDesignMatrix());

    if (globalopts.voxelwise_ev_numbers.size()>0 && globalopts.voxelwiseEvFilenames.size()>0)
      paradigm.loadVoxelwise(globalopts.voxelwise_ev_numbers,globalopts.voxelwiseEvFilenames,mask);
    

    OUT(paradigm.getDesignMatrix().Nrows());
    OUT(paradigm.getDesignMatrix().Ncols());
    OUT(sizeTS);
    OUT(numTS);

    // Setup GLM:
    int numParams = paradigm.getDesignMatrix().Ncols();
    GlimGls glimGls(numTS, sizeTS, numParams);

    // Residuals container:
    Matrix residuals(sizeTS, numTS);

    // Setup autocorrelation estimator:
    AutoCorrEstimator acEst(residuals);

    acEst.mask=mask;

    if(!globalopts.noest)
      {
	cout << "Calculating residuals..." << endl; 
	for(int i = 1; i <= numTS; i++)
	  {						    
            glimGls.setData(datam.Column(i), paradigm.getDesignMatrix(i), i);
	    residuals.Column(i)=glimGls.getResiduals();
	  }
	cout << "Completed" << endl; 
	
	cout << "Estimating residual autocorrelation..." << endl; 
		
	if(globalopts.fitAutoRegressiveModel)
	  {
	    volume4D<float> beta;
	    beta.setmatrix(acEst.fitAutoRegressiveModel(),mask);
	    copybasicproperties(reference,beta);
	    beta.setDisplayMaximumMinimum(beta.max(),beta.min());
	    save_volume4D(beta,LogSingleton::getInstance().getDir() + "/betas");
	  }
	else if(globalopts.tukey)
	  {    
	    if(globalopts.tukeysize == 0)
	      globalopts.tukeysize = (int)(2*sqrt(sizeTS))/2;

	    acEst.calcRaw();

	    if(globalopts.smoothACEst)
		acEst.spatiallySmooth(logger.getDir() + "/" + globalopts.epifname, epivol, globalopts.ms, globalopts.epifname, globalopts.epith, reference[0], globalopts.tukeysize);	

		
	    acEst.tukey(globalopts.tukeysize);
	  }
	else if(globalopts.multitaper)
	  {
	    acEst.calcRaw();
	    acEst.multitaper(int(globalopts.multitapersize));
	  }
	else if(globalopts.pava)
	  {
	    acEst.calcRaw();

	    if(globalopts.smoothACEst)
		acEst.spatiallySmooth(logger.getDir() + "/" + globalopts.epifname, epivol, globalopts.ms, globalopts.epifname, globalopts.epith, reference[0]);
	    
	    acEst.pava();
	  }
	    
      }
    cout << "Completed" << endl; 

    cout << "Prewhitening and Computing PEs..." << endl;
    cout << "Percentage done:" << endl;
    int co = 1;

    Matrix mean_prewhitened_dm(paradigm.getDesignMatrix().Nrows(),paradigm.getDesignMatrix().Ncols());
    mean_prewhitened_dm=0;
    for(int i = 1; i <= numTS; i++)
    {	
      Matrix effectiveDesign(paradigm.getDesignMatrix(i));
      if ( (100.0*i)/numTS > co )
        cout << co++ << "," << flush;	   
      if(!globalopts.noest) {
	acEst.setDesignMatrix(effectiveDesign);
	// Use autocorr estimate to prewhiten data and design:
	ColumnVector xw;
	acEst.preWhiten(datam.Column(i), xw, i, effectiveDesign);   
        datam.Column(i)=xw;
      }
      glimGls.setData(datam.Column(i), effectiveDesign, i);
      residuals.Column(i)=glimGls.getResiduals();
      if(globalopts.output_pwdata || globalopts.verbose)
        mean_prewhitened_dm+=effectiveDesign;	
    }

    if(globalopts.output_pwdata || globalopts.verbose) 
      mean_prewhitened_dm/=numTS;
     
    cerr << "Completed" << endl << "Saving results... " << endl;

    input_data.setmatrix(residuals,mask);
    input_data.setDisplayMaximumMinimum(input_data.max(),input_data.min());
    save_volume4D(input_data,logger.getDir() + "/res4d");

    if(globalopts.output_pwdata || globalopts.verbose)
      {
	// Write out whitened data
        input_data.setmatrix(datam,mask);
	input_data.setDisplayMaximumMinimum(input_data.max(),input_data.min());
        save_volume4D(input_data,logger.getDir() + "/prewhitened_data");
	// Write out whitened design matrix
	write_vest(logger.appendDir("mean_prewhitened_dm.mat"), mean_prewhitened_dm);
		
      }

    // Write out threshac:
    Matrix& threshacm = acEst.getEstimates();
    int cutoff = sizeTS/2;
    if(globalopts.tukey) cutoff = globalopts.tukeysize;
    threshacm = threshacm.Rows(1,MISCMATHS::Max(1,cutoff)); 

    input_data.setmatrix(threshacm,mask);
    input_data.settdim(reference.tdim()); //Possibly just set to a constant 1?
    input_data.set_intent(NIFTI_INTENT_ESTIMATE,0,0,0);
    input_data.setDisplayMaximumMinimum(input_data.max(),input_data.min());
    save_volume4D(input_data,logger.getDir() + "/threshac1");

    threshacm.CleanUp();

    // save gls results:
    glimGls.Save(mask,reference.tdim());
    glimGls.CleanUp();

    cerr << "Completed" << endl;
  }  
  catch(Exception p_excp) 
  {
      cerr << p_excp.what() << endl;
      freeparser(argc, argv);
      return 1;
  }
  catch(...)
  {
      cerr << "Uncaught exception!" << endl;
      freeparser(argc, argv);
      return 1;
  }
  freeparser(argc, argv);
  return 0;
}

