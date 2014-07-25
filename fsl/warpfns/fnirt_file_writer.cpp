// Definitions of class used to write fnirt
// displacement files.
//
// fnirt_file_writer.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2007 University of Oxford 
//

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .save_orig_...
#endif

#include "newimage/newimageall.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "warpfns.h"
#include "fnirt_file_reader.h"
#include "fnirt_file_writer.h"

using namespace std;
using namespace NEWMAT;
using namespace BASISFIELD;
using namespace boost;

namespace NEWIMAGE {

FnirtFileWriter::FnirtFileWriter(const string&                            fname, 
                                 const vector<boost::shared_ptr<basisfield> >&   fields,
                                 Matrix                                   aff)
{
  if (fields.size() != 3 || !fields[0] || !fields[1] || !fields[2]) {
    throw FnirtFileWriterException("FnirtFileWriter: Invalid vector fields");
  }
  common_coef_construction(fname,*(fields[0]),*(fields[1]),*(fields[2]),aff);
}

FnirtFileWriter::FnirtFileWriter(const string&              fname,
				 const volume<float>&       fieldx,                         
			         const volume<float>&       fieldy,                         
			         const volume<float>&       fieldz)
{
  volume<float>  tmp(fieldx.xsize(),fieldx.ysize(),fieldx.zsize());
  tmp.copyproperties(fieldx);
  Matrix         aff = IdentityMatrix(4);

  common_field_construction(fname,tmp,fieldx,fieldy,fieldz,aff);
}

FnirtFileWriter::FnirtFileWriter(const string&                fname,
                                 const volume<float>&         ref,
                                 const volume4D<float>&       vfields,                         
                                 Matrix                       aff)
{
  if (vfields.tsize() != 3) throw FnirtFileWriterException("FnirtFileWriter: Invalid 4D volume");
  common_field_construction(fname,ref,vfields[0],vfields[1],vfields[2],aff);
}

FnirtFileWriter::FnirtFileWriter(const string&                fname,
                                 const volume4D<float>&       vfields)
{
  if (vfields.tsize() != 3) throw FnirtFileWriterException("FnirtFileWriter: Invalid 4D volume");
  volume<float>  tmp(vfields[0].xsize(),vfields[0].ysize(),vfields[0].zsize());
  tmp.copyproperties(vfields[0]);
  Matrix         aff = IdentityMatrix(4);

  common_field_construction(fname,tmp,vfields[0],vfields[1],vfields[2],aff);
}


////////////////////////////////////////////////////////////////////////////////////////////////
//
// The coefficients are saved in a slightly dodgy format where the fields of the niftii
// header are used to store information that is neccessary for us to reconstruct the 
// displacement fields from the coefficients. E.g. the affine starting guess is stored
// in the sform, the knot-spacings in the "pixdims" and the matrix size of the field/template
// in the offsets of the qform.
// For tis reason the coefficient files need to be read written usinge the _orig_ forms
// of the read/write_volume functions. Otherwise the i/o-functions would interpret these
// fields and potentially left-right swap on read/write.
//
////////////////////////////////////////////////////////////////////////////////////////////////

void FnirtFileWriter::common_coef_construction(const string&            fname, 
                                               const basisfield&        fieldx,
                                               const basisfield&        fieldy,
                                               const basisfield&        fieldz,
                                               const Matrix&            aff)
{
  volume4D<float>       coefs(fieldx.CoefSz_x(),fieldx.CoefSz_y(),fieldx.CoefSz_z(),3);
  vector<float>         ksp(3,1.0);
  try {
    const splinefield& f = dynamic_cast<const splinefield& >(fieldx);
    ksp[0] = float(f.Ksp_x()); ksp[1] = float(f.Ksp_y()); ksp[2] = float(f.Ksp_z());
    if (f.Order() == 2) coefs.set_intent(FSL_QUADRATIC_SPLINE_COEFFICIENTS,f.Vxs_x(),f.Vxs_y(),f.Vxs_z());
    else if (f.Order() == 3) coefs.set_intent(FSL_CUBIC_SPLINE_COEFFICIENTS,f.Vxs_x(),f.Vxs_y(),f.Vxs_z());
  }
  catch (...) {  
    try {
      const dctfield& f = dynamic_cast<const dctfield& >(fieldx);
      coefs.set_intent(FSL_DCT_COEFFICIENTS,f.Vxs_x(),f.Vxs_y(),f.Vxs_z());
      throw FnirtFileWriterException("common_coef_construction: Saving of DCT coefficients not yet implemented");
    }
    catch (...) {
      throw FnirtFileWriterException("common_coef_construction: Unknown field type");
    }
  } 
  coefs.setxdim(ksp[0]); coefs.setydim(ksp[1]); coefs.setzdim(ksp[2]);
  Matrix  qform(4,4); 
  qform = IdentityMatrix(4);
  qform(1,4) = fieldx.FieldSz_x(); 
  qform(2,4) = fieldx.FieldSz_y(); 
  qform(3,4) = fieldx.FieldSz_z();
  coefs.set_qform(NIFTI_XFORM_SCANNER_ANAT,qform);
  coefs.set_sform(NIFTI_XFORM_SCANNER_ANAT,aff);
  vector<boost::shared_ptr<ColumnVector> > coefp(3);
  coefp[0]=fieldx.GetCoef(); coefp[1]=fieldy.GetCoef(); coefp[2]=fieldz.GetCoef();  
  for (unsigned int v=0; v<3; v++) {
    ColumnVector  c = *(coefp[v]);
    for (unsigned int k=0, vindx=0; k<fieldx.CoefSz_z(); k++) {
      for (unsigned int j=0; j<fieldx.CoefSz_y(); j++) {
        for (unsigned int i=0; i<fieldx.CoefSz_x(); i++) {
          coefs(i,j,k,v) = c.element(vindx++);
	}
      }
    }
  }

  save_orig_volume4D(coefs,fname);
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
// This saves out the deformation/displacement fields. The affine "starting-guess" will 
// be incorporated into these warps. That means that the header is free to store its "proper"
// information. The fields are in mm and should be added to the "scaled flirt coordinates". 
//
////////////////////////////////////////////////////////////////////////////////////////////////

void FnirtFileWriter::common_field_construction(const string&            fname,
                                                const volume<float>&     ref,
                                                const volume<float>&     fieldx,
                                                const volume<float>&     fieldy,
                                                const volume<float>&     fieldz, 
                                                const Matrix&            aff)
{
  volume4D<float>   fields(ref.xsize(),ref.ysize(),ref.zsize(),3);
  fields.copyproperties(ref); 

  Matrix M;
  bool   add_affine = false;
  if (add_affine = ((aff-IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6)) { // Add affine part to fields
    M = (aff.i() - IdentityMatrix(4))*ref.sampling_mat();
  }

  if (samesize(ref,fieldx,true)) { // If ref is same size as the original field
    fields[0] = fieldx; fields[1] = fieldy; fields[2] = fieldz;
    fields.copyproperties(ref);    // Put qform/sform and stuff back.
    if (add_affine) {
      ColumnVector xv(4), xo(4);
      int zs = ref.zsize(), ys = ref.ysize(), xs = ref.xsize();
      xv(4) = 1.0;
      for (int z=0; z<zs; z++) {
        xv(3) = double(z);
        for (int y=0; y<ys; y++) {
          xv(2) = double(y);
          for (int x=0; x<xs; x++) {
            xv(1) = double(x);
            xo = M*xv;
            fields(x,y,z,0) += xo(1);
            fields(x,y,z,1) += xo(2);
            fields(x,y,z,2) += xo(3);
	  }
        }
      }
    } 
  }
  else {
    fieldx.setextrapolationmethod(extraslice);
    fieldy.setextrapolationmethod(extraslice);
    fieldz.setextrapolationmethod(extraslice);
    Matrix R2F = fieldx.sampling_mat().i() * ref.sampling_mat();
    ColumnVector xv(4), xo(4), xr(4);
    int zs = ref.zsize(), ys = ref.ysize(), xs = ref.xsize();
    xv(4) = 1.0;
    for (int z=0; z<zs; z++) {
      xv(3) = double(z);
      for (int y=0; y<ys; y++) {
        xv(2) = double(y);
        for (int x=0; x<xs; x++) {
          xv(1) = double(x);
          xr = R2F*xv;
          fields(x,y,z,0) = fieldx.interpolate(xr(1),xr(2),xr(3));
          fields(x,y,z,1) = fieldy.interpolate(xr(1),xr(2),xr(3));
          fields(x,y,z,2) = fieldz.interpolate(xr(1),xr(2),xr(3));
          if (add_affine) {
            xo = M*xv;
            fields(x,y,z,0) += xo(1);
            fields(x,y,z,1) += xo(2);
            fields(x,y,z,2) += xo(3);
	  }
        }
      }
    }
  }

  fields.set_intent(FSL_FNIRT_DISPLACEMENT_FIELD,fields.intent_param(0),fields.intent_param(1),fields.intent_param(2));
  fields.setDisplayMaximum(0.0);
  fields.setDisplayMinimum(0.0);

  // Save resulting field
  save_volume4D(fields,fname);
}

/*
void FnirtFileWriter::common_field_construction(const string&            fname,
                                                const volume<float>&     ref,
                                                const volume<float>&     fieldx,
                                                const volume<float>&     fieldy,
                                                const volume<float>&     fieldz, 
                                                const Matrix&            aff)
{
  volume4D<float>   fields(fieldx.xsize(),fieldx.ysize(),fieldx.zsize(),3);
  fields[0] = fieldx; fields[1] = fieldy; fields[2] = fieldz;
  fields.copyproperties(ref); 

  if ((aff-IdentityMatrix(4)).MaximumAbsoluteValue() > 1e-6) {
    // Add affine part to fields
    Matrix M = (aff.i() - IdentityMatrix(4))*ref.sampling_mat();
    ColumnVector xv(4), xo(4);
    int zs = ref.zsize(), ys = ref.ysize(), xs = ref.xsize();
    xv(4) = 1.0;
    for (int z=0; z<zs; z++) {
      xv(3) = double(z);
      for (int y=0; y<ys; y++) {
        xv(2) = double(y);
        for (int x=0; x<xs; x++) {
          xv(1) = double(x);
          xo = M*xv;
          fields(x,y,z,0) += xo(1);
          fields(x,y,z,1) += xo(2);
          fields(x,y,z,2) += xo(3);
        }
      }
    }
  }
  fields.set_intent(FSL_FNIRT_DISPLACEMENT_FIELD,fields.intent_param(0),fields.intent_param(1),fields.intent_param(2));

  // Save resulting field
  save_volume4D(fields,fname);
}
*/

} // End of namespace NEWIMAGE
