//     fslhd.cc - show image header
//     Steve Smith, Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2005 University of Oxford  
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

#include "newimage/newimageall.h"
#include <iomanip>
#include <iostream>
#include "parser.h"

using namespace NEWIMAGE;

namespace fslhd {
void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: fslhd [-x] <input>" << endl;
  cout << "       -x : instead print an XML-style NIFTI header" << endl;
}

void ShowNifti(char *fileName, FSLIO* fslio)
{
    ofstream osaida;
	char filename[500];
    sprintf(filename, "%s.txt", FslMakeBaseName(fileName));
  
    osaida.open(filename);
    if (fslio == NULL) {
        cerr << "ERROR: Could not open file" << endl;
	return;
    }
    if (fslio->niftiptr!=NULL) {
        osaida << nifti_image_to_ascii(fslio->niftiptr) << endl;
	return;
    }
    if (fslio->mincptr!=NULL) {
        osaida << "ERROR: Minc is not currently supported" << endl;
	return;
    }
    return;
}

void ShowHdr(char *fileName, FSLIO* fslio)
{
  int i, ft, isanalyze=0;
  struct dsr *hdr;
  char filename[500];
  mat44 mat;
  int icode, jcode, kcode;
  ofstream osaida;
  sprintf(filename, "%s.txt", FslMakeBaseName(fileName));
  
  osaida.open(filename);

  if (fslio == NULL) 
  {
    cerr << "ERROR: Could not open file" << endl;;
    return;
  }

  ft = FslGetFileType(fslio);
  if (FslBaseFileType(ft)==FSL_TYPE_MINC) 
  {
    cerr << "ERROR: Minc is not currently supported" << endl;
    return;
  }

  if (fslio->niftiptr == NULL) 
  {
    cerr << "ERROR: Not an Analyze or Nifti file" << endl;
    return;
  }
  //-------------------- ANALYZE CASE ----------------------- 
  if (FslBaseFileType(ft)==FSL_TYPE_ANALYZE) 
  { 
    isanalyze=1; 
    //load raw hdr structure   
    hdr = (struct dsr *)calloc(1,sizeof(struct dsr));
    FslReadRawHeader(hdr,fslio->niftiptr->fname);
    if (fslio->niftiptr->byteorder != nifti_short_order()) 
    {
      osaida << "Byte swapping" << endl;
      AvwSwapHeader(hdr);
    }
    osaida << "filename       " << fileName << endl << endl;
    
    //Header Key
    osaida << "sizeof_hdr     " << hdr->hk.sizeof_hdr<< endl;
    osaida << "data_type      " << hdr->hk.data_type << endl;
    osaida << "db_name        " << hdr->hk.db_name << endl;
    osaida << "extents        " << hdr->hk.extents << endl;
    osaida << "session_error  " << hdr->hk.session_error << endl;
    osaida << "regular        " << hdr->hk.regular << endl;
    osaida << "hkey_un0       " << hdr->hk.hkey_un0 << endl;
    //Image Dimension 
    for(i=0;i<8;i++) osaida << "dim" << i << "           " << hdr->dime.dim[i] << endl;
    osaida << "vox_units      " << hdr->dime.vox_units << endl;
    osaida << "cal_units      " << hdr->dime.cal_units << endl;
    osaida << "unused1        " << hdr->dime.unused1 << endl;
    osaida << "datatype       " << hdr->dime.datatype << endl;
    osaida << "bitpix         " << hdr->dime.bitpix << endl;
    osaida.setf(ios::fixed);  //need osaida.setf(ios::fixed) instead of << fixed in stream for tru64 comp
    for(i=0;i<8;i++) osaida << "pixdim" << i << "        " << setprecision(10) << hdr->dime.pixdim[i] << endl;       
    osaida.precision(4);                                   
    osaida << "vox_offset     " << setw(6) << hdr->dime.vox_offset << endl;
    osaida << "funused1       " << setw(6) << hdr->dime.funused1 << endl;
    osaida << "funused2       " << setw(6) << hdr->dime.funused2 << endl;
    osaida << "funused3       " << setw(6) << hdr->dime.funused3 << endl;
    osaida << "cal_max        " << setw(6) << hdr->dime.cal_max << endl;
    osaida << "cal_min        " << setw(6) << hdr->dime.cal_min << endl;
    osaida << "compressed     " << hdr->dime.compressed << endl;
    osaida << "verified       " << hdr->dime.verified << endl;
    osaida << "glmax          " << hdr->dime.glmax << endl;
    osaida << "glmin          " << hdr->dime.glmin << endl ;
    //Data History 
    osaida << "descrip        " <<  hdr->hist.descrip << endl;
    osaida << "aux_file       " <<  hdr->hist.aux_file << endl;
    osaida << "orient         " << (int)hdr->hist.orient << endl; //need cast else blank
    osaida << "originator     " << hdr->hist.originator << endl;
    /*{
      short blah[5];
      memcpy(blah,hdr->hist.originator,5*sizeof(short));
      osaida << "origin1        " << blah[0] << endl;
      osaida << "origin2        " <<blah[1] << endl;
      osaida << "origin3        " << blah[2] << endl;
      }*/
    osaida << "origin1        " << (short)hdr->hist.originator[0] << endl; //These lines don't work on
    osaida << "origin2        " << (short)hdr->hist.originator[2] << endl; //all platforms... but WHICH 
    osaida << "origin3        " << (short)hdr->hist.originator[4] << endl; //ones - that is the question
    osaida << "generated      " << hdr->hist.generated << endl;
    osaida << "scannum        " << hdr->hist.scannum << endl;
    osaida << "patient_id     " <<  hdr->hist.patient_id << endl;
    osaida << "exp_date       " << hdr->hist.exp_date << endl;
    osaida << "exp_time       " << hdr->hist.exp_time << endl;
    osaida << "hist_un0       " << hdr->hist.hist_un0 <<endl;
    osaida << "views          " << hdr->hist.views << endl;
    osaida << "vols_added     " << hdr->hist.vols_added << endl;
    osaida << "start_field    " << hdr->hist.start_field << endl;
    osaida << "field_skip     " << hdr->hist.field_skip << endl;
    osaida << "omax           " << hdr->hist.omax << endl;
    osaida << "omin           " << hdr->hist.omin << endl;
    osaida << "smin           " << hdr->hist.smax << endl;
    osaida << "smin           " << hdr->hist.smin << endl;
    osaida << "file_type      " << FslFileTypeString(0) << endl;
    osaida << "file_code      0" << endl;
    return;
  }
    /* -------------------- NIFTI CASE ----------------------- */
  if (fslio->niftiptr->byteorder != nifti_short_order()) 
  { 
    osaida << "Byte swapping" << endl;
  }

  osaida << "filename       " <<  fslio->niftiptr->fname << endl << endl;
  osaida << "sizeof_hdr     " <<  "348" << endl;
  osaida << "data_type      " <<  nifti_datatype_string(fslio->niftiptr->datatype) << endl;
  for(i=0;i<8;i++)  osaida << "dim" << i << "           " <<  fslio->niftiptr->dim[i] << endl;
  osaida << "vox_units      " <<  nifti_units_string(fslio->niftiptr->xyz_units) << endl;
  osaida << "time_units     " <<  nifti_units_string(fslio->niftiptr->time_units) << endl;
  osaida << "datatype       " <<  fslio->niftiptr->datatype << endl;
  osaida << "nbyper         " <<  fslio->niftiptr->nbyper << endl;
  osaida << "bitpix         " <<  fslio->niftiptr->nbyper * 8 << endl;
  osaida.setf(ios::fixed);  //need osaida.setf(ios::fixed) instead of << fixed in stream for tru64 comp
  for(i=0;i<8;i++) osaida << "pixdim" << i << "        " << setprecision(10) <<  fslio->niftiptr->pixdim[i] << endl;
  osaida << "vox_offset     " <<   fslio->niftiptr->iname_offset << endl;
  osaida.precision(4);
  osaida << "cal_max        " << setw(6) <<  fslio->niftiptr->cal_max << endl;
  osaida << "cal_min        " << setw(6) <<  fslio->niftiptr->cal_min << endl;
  osaida.precision(6);
  osaida << "scl_slope      " <<  fslio->niftiptr->scl_slope << endl;
  osaida << "scl_inter      " <<  fslio->niftiptr->scl_inter << endl;
  osaida << "phase_dim      " <<  fslio->niftiptr->phase_dim << endl;
  osaida << "freq_dim       " <<  fslio->niftiptr->freq_dim  << endl;
  osaida << "slice_dim      " <<  fslio->niftiptr->slice_dim << endl;
  osaida << "slice_name     " <<  nifti_slice_string(fslio->niftiptr->slice_code) << endl;
  osaida << "slice_code     " <<  fslio->niftiptr->slice_code << endl;
  osaida << "slice_start    " <<  fslio->niftiptr->slice_start << endl;
  osaida << "slice_end      " <<  fslio->niftiptr->slice_end << endl;
  osaida << "slice_duration " <<  fslio->niftiptr->slice_duration << endl;
  osaida << "time_offset    " <<  fslio->niftiptr->toffset << endl;
  osaida << "intent         " <<  nifti_intent_string(fslio->niftiptr->intent_code) << endl;
  osaida << "intent_code    " <<  fslio->niftiptr->intent_code << endl;
  osaida << "intent_name    " <<  fslio->niftiptr->intent_name << endl;
  osaida << "intent_p1      " <<  fslio->niftiptr->intent_p1 << endl;
  osaida << "intent_p2      " <<  fslio->niftiptr->intent_p2 << endl;
  osaida << "intent_p3      " <<  fslio->niftiptr->intent_p3 << endl;
  osaida << "qform_name     " <<  nifti_xform_string(fslio->niftiptr->qform_code) << endl;
  osaida << "qform_code     " <<  fslio->niftiptr->qform_code << endl;
  mat = fslio->niftiptr->qto_xyz;
  for(i=1;i<=4;i++) osaida << "qto_xyz:" << i << "      " << mat.m[i-1][0] << "  " << mat.m[i-1][1] << "  " << mat.m[i-1][2] << "  " << mat.m[i-1][3] << endl;
  nifti_mat44_to_orientation(mat,&icode,&jcode,&kcode);
  osaida << "qform_xorient  " << nifti_orientation_string(icode) << endl;
  osaida << "qform_yorient  " << nifti_orientation_string(jcode) << endl;
  osaida << "qform_zorient  " << nifti_orientation_string(kcode) << endl;
  osaida << "sform_name     " << nifti_xform_string(fslio->niftiptr->sform_code) << endl;
  osaida << "sform_code     " << fslio->niftiptr->sform_code << endl;
  mat = fslio->niftiptr->sto_xyz;
  for(i=1;i<=4;i++) osaida << "sto_xyz:" << i << "      " << mat.m[i-1][0] << "  " << mat.m[i-1][1] << "  " << mat.m[i-1][2] << "  " << mat.m[i-1][3] << endl;
  nifti_mat44_to_orientation(mat,&icode,&jcode,&kcode);
  osaida << "sform_xorient  " << nifti_orientation_string(icode) << endl;
  osaida << "sform_yorient  " << nifti_orientation_string(jcode) << endl;
  osaida << "sform_zorient  " << nifti_orientation_string(kcode) << endl;
  osaida << "file_type      " << FslFileTypeString(fslio->niftiptr->nifti_type) << endl;
  osaida << "file_code      " << fslio->niftiptr->nifti_type << endl;
  //Data History 
  osaida << "descrip        " << fslio->niftiptr->descrip << endl;
  osaida << "aux_file       " << fslio->niftiptr->aux_file << endl;
  /*osaida << "orient         %d\n", hdr->hist.orient); */
  return;
}


extern "C" __declspec(dllexport) int _stdcall fslhd(char *CmdLn)
{
  int argc;
  char **argv;
  
  parser(CmdLn, argc, argv);
  if (argc < 2) 
  {
    print_usage(string(argv[0]));
    freeparser(argc, argv);
    return 1; 
  }
  FSLIO* fslio=NULL;
  int argval=1, niftiform=0;
  if (strcmp(argv[1],"-x")==0) 
  {
      niftiform=1;
      argval=2;
  }
  fslio = FslOpen(FslMakeBaseName(argv[argval]),"rb");
  FslClose(fslio);
  if (niftiform==0) ShowHdr(argv[argval], fslio); 
  else ShowNifti(argv[argval],fslio); 
  freeparser(argc, argv);
  return 0;
}


}