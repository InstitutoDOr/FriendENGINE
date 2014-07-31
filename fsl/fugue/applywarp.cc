
/*  applywarp.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

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

#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "warpfns/warpfns.h"
#include "warpfns/fnirt_file_reader.h"
#include "parser.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

const int LargeIma = 512*512*512;    // Used to detect silly super-sampling

using namespace std;
using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace fslapplywarp {
// Does the job
int applywarp();

// Downsamples supersampled output image.
void downsample(const volume<float>&           ivol,
                const vector<unsigned int>&    ss,
		bool                           nn,
                volume<float>&                 ovol);

// Returns most common value in vector
float hist_mode(vector<float>  vec);



////////////////////////////////////////////////////////////////////////////

// COMMAND LINE OPTIONS

string title="applywarp (Version 1.2)\nCopyright(c) 2001, University of Oxford (Mark Jenkinson)";
string examples=string("applywarp -i invol -o outvol -r refvol -w warpvol\n") +
                string("applywarp -i invol -o outvol -r refvol -w coefvol\n");

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> abswarp(string("--abs"), false,
		  string("\ttreat warp field as absolute: x' = w(x)"),
		  false, no_argument);
Option<bool> relwarp(string("--rel"), false,
		  string("\ttreat warp field as relative: x' = x + w(x)"),
		  false, no_argument);
Option<string> interp(string("--interp"), string(""),
		   string("interpolation method {nn,trilinear,sinc,spline}"),
		   false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		       string("filename of input image (to be warped)"),
		       true, requires_argument);
Option<string> refname(string("-r,--ref"), string(""),
		       string("filename for reference image"),
		       true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("filename for output (warped) image"),
		       true, requires_argument);
Option<bool> supersample(string("-s,--super"),false,
			string("intermediary supersampling of output, default is off"),
                        false, no_argument);
Option<string> supersamplelevel(string("--superlevel"),string("2"),
				string("level of intermediary supersampling, a for 'automatic' or integer level. Default = 2"),
				false, requires_argument);
Option<string> datatype(string("-d,--datatype"), string(""),
                        string("Force output data type [char short int float double]."),
                        false, requires_argument);
Option<string> warpname(string("-w,--warp"), string(""),
			string("filename for warp/coefficient (volume)"),
			false, requires_argument);
Option<string> maskname(string("-m,--mask"), string(""),
		       string("filename for mask image (in reference space)"),
		       false, requires_argument);
Option<string> prematname(string("--premat"), string(""),
		       string("filename for pre-transform (affine matrix)"),
		       false, requires_argument);
Option<string> postmatname(string("--postmat"), string(""),
		       string("filename for post-transform (affine matrix)"),
		       false, requires_argument);


int applywarp()
{
  // Check for sillines
  if (relwarp.value() && abswarp.value()) {
    cerr << "Only one of --abs and --rel can be set" << endl;
    exit(EXIT_FAILURE);
  }
  // Assert value for data-type
  short  dtypecode = DT_FLOAT;
  if (datatype.set()) {
    if (datatype.value()==string("char")) dtypecode = DT_UNSIGNED_CHAR;
    else if (datatype.value()==string("short")) dtypecode = DT_SIGNED_SHORT;
    else if (datatype.value()==string("int")) dtypecode = DT_SIGNED_INT;
    else if (datatype.value()==string("float")) dtypecode = DT_FLOAT;
    else if (datatype.value()==string("double")) dtypecode = DT_DOUBLE;
    else {
      cerr << "Unknown data type " << datatype.value() << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  // read in pre/post transforms
  Matrix premat, postmat;
  premat = IdentityMatrix(4);
  postmat = IdentityMatrix(4);
  if (prematname.set()) {
    premat = read_ascii_matrix(prematname.value());
  }
  if (postmatname.set()) {
    postmat = read_ascii_matrix(postmatname.value());
  }

  // read in-images
  volume4D<float> invol;
  read_volume4D(invol,inname.value());
  
  //
  // Get size of output from --ref. 
  //
  volume<float>    refvol;
  read_volume(refvol,refname.value());
  volume4D<float>  outvol;
  for (int i=0; i<invol.tsize(); i++) outvol.addvolume(refvol);

  // Read in/create warps from file
  FnirtFileReader  fnirtfile;
  AbsOrRelWarps    wt = UnknownWarps;
  volume4D<float>  defvol(refvol.xsize(),refvol.ysize(),refvol.zsize(),3);
  Matrix           affmat(4,4);
  if (warpname.set()) {
    if (abswarp.value()) wt = AbsoluteWarps;
    else if (relwarp.value()) wt = RelativeWarps;
    try {
      fnirtfile.Read(warpname.value(),wt,verbose.value());
      defvol = fnirtfile.FieldAsNewimageVolume4D();
      affmat = fnirtfile.AffineMat();
    }
    catch (...) {
      cerr << "An error occured while reading file: " << warpname.value() << endl;
      exit(EXIT_FAILURE);
    }
  }
  else {  // This is intended to let one use applywarp to resample files using flirt-matrices only in a convenient way
    wt = RelativeWarps;
    defvol.setdims(refvol.xdim(),refvol.ydim(),refvol.zdim(),1.0);
    defvol = 0.0;
    affmat = IdentityMatrix(4);
  }

  //
  // Assert and decode supersampling parameters
  //
  vector<unsigned int> ssvec(3,0);
  bool superflag = false;
  if (supersample.value() || supersamplelevel.set()) {
    superflag = true;
    if (supersamplelevel.value() == string("a")) { // If "automatic" supersampling
      ssvec[0] = static_cast<unsigned int>(refvol.xdim() / invol.xdim() + 0.9);
      ssvec[0] = (ssvec[0] > 0) ? ssvec[0] : 1;
      ssvec[1] = static_cast<unsigned int>(refvol.ydim() / invol.ydim() + 0.9);
      ssvec[1] = (ssvec[1] > 0) ? ssvec[1] : 1;
      ssvec[2] = static_cast<unsigned int>(refvol.zdim() / invol.zdim() + 0.9);
      ssvec[2] = (ssvec[2] > 0) ? ssvec[2] : 1;
    }
    else {
      unsigned int ssfac = 0;
      char         skrutt[256];
      if (sscanf(supersamplelevel.value().c_str(),"%1u%s",&ssfac,skrutt) != 1) {
        cerr << "Invalid argument " << supersamplelevel.value() << " to --superlevel parameter" << endl;
        exit(EXIT_FAILURE);
      }
      if (ssfac < 1 || ssfac > 10) { 
        cerr << "Argument to --superlevel parameter must be between 1 and 10, or a (for automatic)" << endl;
        exit(EXIT_FAILURE);
      }
      ssvec.assign(3,ssfac);
    }
  }
  if (verbose.value()) {
    cout << "superflag = " << superflag << endl;
    cout << "ssvec = " << ssvec[0] << "  " << ssvec[1] << "  " << ssvec[2] << endl;
  }

  //
  // Read and verify mask
  //
  volume<float>    mask;
  if (maskname.set()) { 
    read_volume(mask,maskname.value());
    if (!samesize(refvol,mask)) {
      cerr << "--ref and --mask must have same size" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // set interpolation method
  if (interp.value() == "nn" ) {
    invol.setinterpolationmethod(nearestneighbour);
  } 
  else if (!interp.set() || interp.value() == "trilinear") {
    invol.setinterpolationmethod(trilinear);
  } 
  else if (interp.value() == "sinc") {
    invol.setinterpolationmethod(sinc);
  }
  else if (interp.value() == "spline") {
    invol.setinterpolationmethod(spline);
    invol.setsplineorder(3);
  }
  else {
    cerr << "Unknown interpolation type " << interp.value() << endl;
    exit(EXIT_FAILURE);
  }

  boost::shared_ptr<volume<float> >   ssout_ptr;     // Used for (optional) super-sampling
  if (superflag) {
    // Make temporary image volume for use when resampling
    vector<int>   isize(3,0);
    isize[0] = ssvec[0]*outvol.xsize();
    isize[1] = ssvec[1]*outvol.ysize();
    isize[2] = ssvec[2]*outvol.zsize();
    if ((isize[0]*isize[1]*isize[2]) > LargeIma) {
      cerr << "Supersampling renders output image too large" << endl;
      exit(EXIT_FAILURE);
    }
    ssout_ptr = boost::shared_ptr<volume<float> >(new volume<float>(isize[0],isize[1],isize[2]));
    vector<float>  vsize(3,0.0);
    vsize[0] = outvol.xdim() / float(ssvec[0]);
    vsize[1] = outvol.ydim() / float(ssvec[1]);
    vsize[2] = outvol.zdim() / float(ssvec[2]);
    ssout_ptr->setdims(vsize[0],vsize[1],vsize[2]);
    // Correct for half-voxel shift caused by 0,0,0 mm being
    // set at centre of voxel 0,0,0. A bit fiddly because
    // it has to be incorporated into postmat.
    Matrix  M_translate = IdentityMatrix(4);
    M_translate(1,4) = - (outvol.xdim()/float(2) - outvol.xdim()/(float(2)*float(ssvec[0])));
    M_translate(2,4) = - (outvol.ydim()/float(2) - outvol.ydim()/(float(2)*float(ssvec[1])));
    M_translate(3,4) = - (outvol.zdim()/float(2) - outvol.zdim()/(float(2)*float(ssvec[2])));
    postmat = (postmat.i() * M_translate).i();
  }

  for (int t=0; t<invol.tsize(); t++) {
    invol[t].setpadvalue(invol[t].backgroundval());
    invol[t].setextrapolationmethod(extraslice);
    // do the deed
    if (superflag) {
      apply_warp(invol[t],affmat,defvol,postmat,premat,*ssout_ptr);
      if (interp.value() == "nn") downsample(*ssout_ptr,ssvec,true,outvol[t]);
      else downsample(*ssout_ptr,ssvec,false,outvol[t]);
    }
    else {
      apply_warp(invol[t],affmat,defvol,postmat,premat,outvol[t]);
    }
    if (maskname.set()) { outvol[t] *= mask; }
  }

  outvol.setDisplayMaximumMinimum(0,0);
  // save the results
  if (datatype.set()) {
    save_volume4D_dtype(outvol,outname.value(),dtypecode);
  }
  else {
    save_volume4D_dtype(outvol,outname.value(),dtype(inname.value()));
  }
  
  return(EXIT_SUCCESS);
}


extern "C" __declspec(dllexport) int _stdcall applywarp(char *CmdLn)
{
  Tracer tr("main");

  int argc;
  char **argv;
  
  OptionParser options(title, examples);

  parser(CmdLn, argc, argv);
  try {
    options.add(inname);
    options.add(refname);
    options.add(warpname);
    options.add(abswarp);
    options.add(relwarp);
    options.add(outname);
    options.add(datatype);
    options.add(supersample);
    options.add(supersamplelevel);
    options.add(prematname);
    options.add(postmatname);
    options.add(maskname);
    options.add(interp);
    options.add(verbose);
    options.add(help);

    int i=options.parse_command_line(argc, argv);
    if (i < argc) {
      for (; i<argc; i++) {
        cerr << "Unknown input: " << argv[i] << endl;
      }
      freeparser(argc, argv);
      return(EXIT_FAILURE);
    }

    if (help.value() || !options.check_compulsory_arguments(true)) {
      options.usage();
      freeparser(argc, argv);
      exit(EXIT_FAILURE);
    }
  }  
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    freeparser(argc, argv);
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  int r =  applywarp();
  freeparser(argc, argv);
  return(r);
}

void downsample(const volume<float>&          ivol,
                const vector<unsigned int>&   ss,
                bool                          nn,
                volume<float>&                ovol)
{
  ovol = 0.0;
  if (!nn) {
    for (int k=0, kk=0; k<ivol.zsize(); k++) {
      for (int j=0, jj=0; j<ivol.ysize(); j++) {
	for (int i=0, ii=0; i<ivol.xsize(); i++) {
	  ovol(ii,jj,kk) += ivol(i,j,k);
	  if (i%ss[0] == (ss[0]-1)) ii++;
	}
	if (j%ss[1] == (ss[1]-1)) jj++;
      }
      if (k%ss[2] == (ss[2]-1)) kk++;
    }
    ovol /= static_cast<float>(ss[0]*ss[1]*ss[2]);
  }
  else {  // If nearest neighbour we do mode instead of mean
    vector<float>  hvals(ss[0]*ss[1]*ss[2],0.0);
    for (unsigned int kk=0; kk<static_cast<unsigned int>(ovol.zsize()); kk++) {
      for (unsigned int jj=0; jj<static_cast<unsigned int>(ovol.ysize()); jj++) {
	for (unsigned int ii=0; ii<static_cast<unsigned int>(ovol.xsize()); ii++) {
	  unsigned int indx = 0;
          for (unsigned int k=ss[2]*kk; k<ss[2]*kk+ss[2]; k++) {
            for (unsigned int j=ss[1]*jj; j<ss[1]*jj+ss[1]; j++) {
              for (unsigned int i=ss[0]*ii; i<ss[0]*ii+ss[0]; i++) {
                hvals[indx] = ivol(i,j,k);
		indx++;
	      }
	    }
	  }
          ovol(ii,jj,kk) = hist_mode(hvals);
	}
      }
    }
  }
}

float hist_mode(vector<float>  vec)
{
  map<float,unsigned int>            hist;
  map<float,unsigned int>::iterator  pos;
  for (unsigned int i=0; i<vec.size(); i++) hist[vec[i]] += 1;     // Generate histogram
  unsigned int maxcnt=0;
  float        modeval=0.0;
  for (pos=hist.begin(); pos!=hist.end(); ++pos) {                 // Find mode of histogram
    if (pos->second > maxcnt) {maxcnt = pos->second; modeval = pos->first;}
  }
  if (maxcnt==1) { // If there is no mode, use median instead
    unsigned int indx=0;
    for (pos=hist.begin(); pos!=hist.end() && indx<((vec.size()-1)/2); ++pos, indx++) ; // Yes, it is intentional
    modeval = pos->first;
    if (!(vec.size()%2)) { // If even
      ++pos;
      modeval = (modeval + pos->first) / 2.0;      
    }
  }
  return(modeval);
}
}