/*  fsl_regfilt - 

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

namespace fslregfilt { 
// The two strings below specify the title and example usage that is
// printed out as the help or usage message

  string title=string("fsl_regfilt (Version 1.0)")+
		string("\n\n Copyright(c) 2008, University of Oxford (Christian F. Beckmann)\n")+
		string(" Data de-noising by regressing out part of a design matrix\n")+
		string(" using simple OLS regression on 4D images");
  string examples="fsl_regfilt -i <input> -d <design> -f <components> -o <out> [options]";

//Command line Options {
  Option<string> fnin(string("-i,--in"), string(""),
		string("        input file name (4D image)"),
		true, requires_argument);
  Option<string> fnout(string("-o,--out"), string(""),
		string("output file name for the filtered data"),
		true, requires_argument);
  Option<string> fndesign(string("-d,--design"), string(""),
		string("file name of the matrix with time courses (e.g. GLM design or MELODIC mixing matrix)"),
		true, requires_argument);
  Option<string> fnmask(string("-m,--mask"), string(""),
		string("mask image file name"),
		false, requires_argument);
	Option<string> filter(string("-f,--filter"),string(""),
		string("filter out part of the regression model, e.g. -f \"1,2,3\""),
		true, requires_argument);
	Option<bool> verbose(string("-v"),FALSE,
		string("        switch on diagnostic messages"),
		false, no_argument);
	Option<bool> perfvn(string("--vn"),FALSE,
		string("        perfrom variance-normalisation on data"),
		false, no_argument);
	Option<int> help(string("-h,--help"), 0,
		string("display this help text"),
		false,no_argument);
	// Output options	
	Option<string> outdata(string("--out_data"),string(""),
		string("output data"),
		false, requires_argument);
	Option<string> outvnscales(string("--out_vnscales"),string(""),
		string("output scaling factors for variance normalisation"),
		false, requires_argument);
		/*
}
*/
//Globals {
	int voxels = 0;
	Matrix data;
	Matrix design;
	Matrix meanR;
	RowVector vnscales;
	volume<float> mask;
	volume<float> Mean;  /*
}
*/
////////////////////////////////////////////////////////////////////////////

// Local functions
void save4D(Matrix what, string fname){
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

void saveit(Matrix what, string fname){
	if(isimage(what))
		save4D(what,fname);
	else
		write_ascii_matrix(what,fname);
}

int dofilter(){
	if(!isimage(data)){
		cerr << "ERROR: need to specify 4D input to use filtering" << endl;
		return 1;
	}
	if(verbose.value())
		cout << "  Calculating maps " << endl;  
	Matrix unmixMatrix = pinv(design);
  Matrix maps = unmixMatrix * data;

  Matrix noisedes;
  Matrix noisemaps;

  int ctr=0;    
  char *p;
  char t[1024];
  const char *discard = ", [];{(})abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*_-=+|\':><./?";
  
  strcpy(t, filter.value().c_str());
  p=strtok(t,discard);
  ctr = atoi(p);
  if(ctr>0 && ctr<=design.Ncols()){
    noisedes = design.Column(ctr);
    noisemaps  = maps.Row(ctr).t();    
  }

  do{
    p=strtok(NULL,discard);
    if(p){
			ctr = atoi(p);
      if(ctr>0 && ctr<=design.Ncols()){
  			noisedes |= design.Column(ctr);
  			noisemaps  |= maps.Row(ctr).t();
			}
    }
  }while(p);
  Matrix newData;
	if(verbose.value())
		cout << "  Calculating filtered data " << endl;
 	newData = data - noisedes * noisemaps.t();
  newData = newData + ones(newData.Nrows(),1)*meanR;
  
	save4D(newData,fnout.value());
  return 0;	
}

int setup(){
	if(fsl_imageexists(fnin.value())){//read data
		//input is 3D/4D vol
		volume4D<float> tmpdata;
		read_volume4D(tmpdata,fnin.value());
		
		// create mask
		if(fnmask.value()>""){
			read_volume(mask,fnmask.value());
			if(!samesize(tmpdata[0],mask)){
				cerr << "ERROR: Mask image does not match input image" << endl;
				return 1;
			};
		}else{
			if(verbose.value())
				cout << "  Creating mask image "  << endl;
			Mean = meanvol(tmpdata);
			float Mmin, Mmax;
			Mmin = Mean.min(); Mmax = Mean.max();
			mask = binarise(Mean,float(Mmin + 0.01* (Mmax-Mmin)),Mmax);
		}
		
		data = tmpdata.matrix(mask);
		voxels = data.Ncols();
		if(verbose.value())
			cout << " Data matrix size : " << data.Nrows() << " x " << voxels << endl;
		
	}else{
		cerr << "ERROR: cannot read input image " << fnin.value()<<endl;
		return 1;
	}

	design = read_ascii_matrix(fndesign.value());

	meanR=mean(data,1);
	data = remmean(data,1);
	design = remmean(design,1);
	if(perfvn.value())
		vnscales = Melodic::varnorm(data);
	return 0;	
}

void write_res(){	
	if(outdata.value()>"")
		saveit(data,outdata.value());
	if(outvnscales.value()>"")
		saveit(vnscales,outvnscales.value());
}

int do_work(int argc, char* argv[]) {
  if(setup())
		exit(1);

	if(dofilter())
		exit(1);	
	write_res();
	return 0;
}

////////////////////////////////////////////////////////////////////////////

extern "C" __declspec(dllexport) int _stdcall fsl_regfilt(char *CmdLn)
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
			options.add(fnmask);
			options.add(filter);
			options.add(perfvn);
			options.add(verbose);
			options.add(help);
			options.add(outdata);
			options.add(outvnscales);
	    options.parse_command_line(argc, argv);

	    // line below stops the program if the help was requested or 
	    //  a compulsory option was not set
	    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
				options.usage();
				r =EXIT_FAILURE;
	    }else{
	  		// Call the local functions
	  	  	    r = do_work(argc,argv);
			}
		}catch(X_OptionError& e) {
			options.usage();
	  	cerr << endl << e.what() << endl;
	    r = EXIT_FAILURE;
	  }catch(std::exception &e) {
	    cerr << e.what() << endl;
	  } 
      freeparser(argc, argv);
      return r;

}

}