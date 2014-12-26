/*  fsl_glm - 

    Christian F. Beckmann, FMRIB Image Analysis Group

    Copyright (C) 2006-2008 University of Oxford  */

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

#include "libvis/miscplot.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "utils/options.h"
#include <vector>
#include "newimage/newimageall.h"
#include "melhlprfns.h"
#include "parser.h"

using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace Utilities;
using namespace std;

namespace fslglm
{
// The two strings below specify the title and example usage that is
// printed out as the help or usage message

  string title=string("fsl_glm (Version 1.1)")+
		string("\nCopyright(c) 2004-2009, University of Oxford (Christian F. Beckmann)\n")+
		string(" \n Simple GLM usign ordinary least-squares (OLS) regression on\n")+
		string(" time courses and/or 3D/4D imges against time courses \n")+
		string(" or 3D/4D images");
  string examples="fsl_glm -i <input> -d <design> -o <output> [options]";

//Command line Options {
  Option<string> fnin(string("-i,--in"), string(""),
		string("        input file name (text matrix or 3D/4D image file)"),
		true, requires_argument);
  Option<string> fnout(string("-o,--out"), string(""),
		string("output file name for GLM parameter estimates (GLM betas)"),
		false, requires_argument);
  Option<string> fndesign(string("-d,--design"), string(""),
		string("file name of the GLM design matrix (time courses or spatial maps)"),
		false, requires_argument);
  Option<string> fnmask(string("-m,--mask"), string(""),
		string("mask image file name if input is image"),
		false, requires_argument);
  Option<string> fncontrasts(string("-c,--contrasts"), string(""),
		string("matrix of t-statistics contrasts"),
		false, requires_argument);
  Option<string> fnftest(string("-f,--ftests"), string(""),
		string("matrix of F-tests on contrasts"),
		false, requires_argument,false);
	Option<int> dofset(string("--dof"),0,
		string("        set degrees-of-freedom explicitly"),
		false, requires_argument);
	Option<bool> normdes(string("--des_norm"),FALSE,
		string("switch on normalisation of the design matrix columns to unit std. deviation"),
		false, no_argument);
	Option<bool> normdat(string("--dat_norm"),FALSE,
		string("switch on normalisation of the data time series to unit std. deviation"),
		false, no_argument);
	Option<bool> perfvn(string("--vn"),FALSE,
		string("        perform MELODIC variance-normalisation on data"),
		false, no_argument);
	Option<bool> perf_demean(string("--demean"),FALSE,
		string("switch on de-meaning of design and data"),
		false, no_argument);
	Option<int> help(string("-h,--help"), 0,
		string("display this help text"),
		false,no_argument);
	Option<bool> debug(string("--debug"), FALSE,
		string("display debug information"),
		false,no_argument,false);
	// Output options	
	Option<string> outcope(string("--out_cope"),string(""),
		string("output file name for COPEs (either as text file or image)"),
		false, requires_argument);
	Option<string> outz(string("--out_z"),string(""),
		string("        output file name for Z-stats (either as text file or image)"),
		false, requires_argument);
	Option<string> outt(string("--out_t"),string(""),
		string("        output file name for t-stats (either as text file or image)"),
		false, requires_argument);
	Option<string> outp(string("--out_p"),string(""),
		string("        output file name for p-values of Z-stats (either as text file or image)"),
		false, requires_argument);
	Option<string> outf(string("--out_f"),string(""),
		string("        output file name for F-value of full model fit"),
		false, requires_argument);
	Option<string> outpf(string("--out_pf"),string(""),
		string("output file name for p-value for full model fit"),
		false, requires_argument);
	Option<string> outres(string("--out_res"),string(""),
		string("output file name for residuals"),
		false, requires_argument);
	Option<string> outvarcb(string("--out_varcb"),string(""),
		string("output file name for variance of COPEs"),
		false, requires_argument);
	Option<string> outsigsq(string("--out_sigsq"),string(""),
		string("output file name for residual noise variance sigma-square"),
		false, requires_argument);
	Option<string> outdata(string("--out_data"),string(""),
		string("output file name for pre-processed data"),
		false, requires_argument);
	Option<string> outvnscales(string("--out_vnscales"),string(""),
		string("output file name for scaling factors for variance normalisation"),
		false, requires_argument);
		/*
}
*/
	int voxels = 0;

////////////////////////////////////////////////////////////////////////////

// Local functions
void save4D(Matrix what, string fname, Matrix &data, volume<float> &mask){
		if(what.Ncols()==data.Ncols()||what.Nrows()==data.Nrows()){
			volume4D<float> tempVol;
			if(what.Nrows()>what.Ncols())
				tempVol.setmatrix(what.t(),mask);
			else
				tempVol.setmatrix(what,mask);
			save_volume4D(tempVol,fname);
		}
}

bool isimage(Matrix what){
	if((voxels > 0)&&(what.Ncols()==voxels || what.Nrows()==voxels))
		return TRUE;
	else
		return FALSE;
}

void saveit(Matrix what, string fname, Matrix &data, volume<float> &mask){
	if(isimage(what))
		save4D(what,fname, data, mask);
	else if(fsl_imageexists(fndesign.value()))
		write_ascii_matrix(what.t(),fname);
	else
		write_ascii_matrix(what,fname);
}

int setup(Matrix &data, volume<float> &mask, Matrix &design, Matrix &contrasts, Matrix &fcontrasts, Matrix &meanR, RowVector &vnscales)
{
	if(fsl_imageexists(fnin.value())){//read data
		//input is 3D/4D vol
		volume4D<float> tmpdata;
		read_volume4D(tmpdata,fnin.value());
		
		// create mask
		if(fnmask.value()>""){
			if(debug.value())
				cout << "Reading mask file " << fnmask.value() << endl;
			read_volume(mask,fnmask.value());
			if(!samesize(tmpdata[0],mask)){
				cerr << "ERROR: Mask image does not match input image" << endl;
				return 1;
			};
		}else{
			if(debug.value())
				cout << "Creating mask image" << endl;
			mask=tmpdata[0]*0.0+1.0;
			data=tmpdata.matrix(mask);
			Melodic::update_mask(mask,data);
		}	
		data = tmpdata.matrix(mask);
		voxels = data.Ncols();
			
		if(perfvn.value()){
			if(debug.value())
				cout << "Perform MELODIC variance normalisation (and demeaning)" << endl;			
			data = remmean(data,1);
			vnscales = Melodic::varnorm(data);		
		}
	}
	else
		data = read_ascii_matrix(fnin.value());	

	if(fsl_imageexists(fndesign.value())){//read design
		if(debug.value())
			cout << "Reading design file "<< fndesign.value()<< endl;
		volume4D<float> tmpdata;
		read_volume4D(tmpdata,fndesign.value());
		if(!samesize(tmpdata[0],mask)){
			cerr << "ERROR: GLM design does not match input image in size" << endl;
			return 1;
		}
		if(debug.value())
			cout << "Transposing data" << endl;
		design = tmpdata.matrix(mask).t();
		data = data.t();
	}else{
		design = read_ascii_matrix(fndesign.value());
	}

	if(perf_demean.value()){
		if(debug.value())
			cout << "De-meaning the data matrix" << endl;
		data = remmean(data,1);
	}
	if(normdat.value()){
		if(debug.value())
			cout << "Normalising data matrix to unit std-deviation" << endl;
		data =  SP(data,ones(data.Nrows(),1)*pow(stdev(data,1),-1));
	}

	meanR=mean(data,1);
	if(perf_demean.value()){
		if(debug.value())
			cout << "De-meaning design matrix" << endl;	
		design = remmean(design,1);
	}
	if(normdes.value()){
		if(debug.value())
			cout << "Normalising design matrix to unit std-deviation" << endl;
		design =  SP(design,ones(design.Nrows(),1)*pow(stdev(design,1),-1));
	}
	if(fncontrasts.value()>""){//read contrast		
		contrasts = read_ascii_matrix(fncontrasts.value());
		if(!(contrasts.Ncols()==design.Ncols())){
			cerr << "ERROR: contrast matrix GLM design does not match GLM design" << endl;
			return 1;
		}
	}else{
		contrasts = IdentityMatrix(design.Ncols());
		contrasts &= -1.0 * contrasts;
	}
	return 0;	
}

void write_res(Melodic::basicGLM &glm, Matrix &data, RowVector &vnscales, volume<float> &mask){
	if(fnout.value()>"")
		saveit(glm.get_beta(),fnout.value(), data, mask);
	if(outcope.value()>"")
		saveit(glm.get_cbeta(), outcope.value(), data, mask);
	if(outz.value()>"")
		saveit(glm.get_z(), outz.value(), data, mask);
	if(outt.value()>"")
		saveit(glm.get_t(), outt.value(), data, mask);
	if(outp.value()>"")
		saveit(glm.get_p(), outp.value(), data, mask);
	if(outf.value()>"")
		saveit(glm.get_f_fmf(), outf.value(), data, mask);
	if(outpf.value()>"")
		saveit(glm.get_pf_fmf(), outpf.value(), data, mask);
	if(outres.value()>"")
		saveit(glm.get_residu(), outres.value(), data, mask);
	if(outvarcb.value()>"")
		saveit(glm.get_varcb(), outvarcb.value(), data, mask);
	if(outsigsq.value()>"")
		saveit(glm.get_sigsq(), outsigsq.value(), data, mask);
	if(outdata.value()>"")
		saveit(data, outdata.value(), data, mask);
	if(outvnscales.value()>"")
		saveit(vnscales, outvnscales.value(), data, mask);
}

int do_work(int argc, char* argv[]) {
	//Globals 
	Melodic::basicGLM glm;
	Matrix data;
	Matrix design;
	Matrix contrasts;
	Matrix fcontrasts;
	Matrix meanR;
	RowVector vnscales;
	volume<float> mask;

	if (setup(data, mask, design, contrasts, fcontrasts, meanR, vnscales))
		return(1);

	glm.olsfit(data,design,contrasts,dofset.value());
	write_res(glm, data, vnscales, mask);
	return 0;
}

////////////////////////////////////////////////////////////////////////////

extern "C" __declspec(dllexport) int _stdcall fsl_glm(char *CmdLn)
{
      int r;
	  int argc;
      char **argv;
  
      parser(CmdLn, argc, argv);
      Tracer tr("main");
	  OptionParser options(title, examples);
	  try{
	    // must include all wanted options here (the order determines how
	    //  the help message is printed)
			options.add(fnin);
			options.add(fnout);
			options.add(fndesign);
			options.add(fncontrasts);
			options.add(fnmask);
			options.add(fnftest);
			options.add(dofset);
			options.add(normdes);
			options.add(normdat);
			options.add(perfvn);
			options.add(perf_demean);
			options.add(help);
			options.add(debug);
			options.add(outcope);
			options.add(outz);
			options.add(outt);
			options.add(outp);
			options.add(outf);
			options.add(outpf);
			options.add(outres);
			options.add(outvarcb);
			options.add(outsigsq);
			options.add(outdata);
			options.add(outvnscales);
	    options.parse_command_line(argc, argv);

	    // line below stops the program if the help was requested or 
	    //  a compulsory option was not set
	    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
				options.usage();
				r=EXIT_FAILURE;
	    }else{
	  		// Call the local functions
	  		r= do_work(argc,argv);
			}
		}catch(X_OptionError& e) {
			options.usage();
	  	cerr << endl << e.what() << endl;
	    r=EXIT_FAILURE;
	  }catch(std::exception &e) {
	    cerr << e.what() << endl;
	  } 
      freeparser(argc, argv);
      return r;
	}

}