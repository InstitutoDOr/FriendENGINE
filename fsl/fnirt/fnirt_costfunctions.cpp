// Definitions of cost-function clases used by FNIRT
//
// fnirt_costfunctions.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2007 University of Oxford 
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .set_sform etc
#endif
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "warpfns/warpfns.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "warpfns/fnirt_file_writer.h"
#include "matching_points.h"
#include "fnirt_costfunctions.h"


#include "fnirtfns.h"  // For debugging


#ifndef SQR
#define SQR(A)  (A)*(A)
#endif

using namespace FNIRT;

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Member functions for class fnirt_CF
//
///////////////////////////////////////////////////////////////////////////////////////////////

fnirt_CF::fnirt_CF(const NEWIMAGE::volume<float>&                                  pvref, 
                   const NEWIMAGE::volume<float>&                                  pvobj,
                   const NEWMAT::Matrix&                                           paff,
                   std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&        pdf_field,
                   boost::shared_ptr<IntensityMapper>                              im)
: vref(new NEWIMAGE::volume<float>(pvref)), vobj(new NEWIMAGE::volume<float>(pvobj)), aff(paff), df_field(3), imap(im)
{
  common_construction(pdf_field);
}

fnirt_CF::fnirt_CF(const boost::shared_ptr<NEWIMAGE::volume<float> >               pvref, 
                   const boost::shared_ptr<NEWIMAGE::volume<float> >               pvobj,
                   const NEWMAT::Matrix&                                           paff,
                   std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&        pdf_field,
                   boost::shared_ptr<IntensityMapper>                              im)
: vref(pvref), vobj(pvobj), aff(paff), df_field(3), imap(im)
{
  common_construction(pdf_field);
}

void fnirt_CF::common_construction(std::vector<boost::shared_ptr<BASISFIELD::basisfield> >& pdf_field)
{
  // Do some checks to see all is kosher

  if (pdf_field.size() != 3) {throw FnirtException("fnirt_CF::fnirt_CF: 3 fields must be input");}
  if (OrigRefSz_x() != pdf_field[0]->FieldSz_x() || OrigRefSz_y() != pdf_field[0]->FieldSz_y() || OrigRefSz_z() != pdf_field[0]->FieldSz_z()) { 
    throw FnirtException("fnirt_CF::fnirt_CF: Mismatch between reference volume and basis field");
  }

  if (!imap) {
    imap = boost::shared_ptr<IntensityMapper>(new IntensityMapper);  // NO_SCALING if no imap was specified
  }

  for (int i=0; i<3; i++) {df_field[i] = pdf_field[i];}

  // Since no smoothing or subsampling has yet to be applied to the ref image we point
  // the smoothed and subsampled version straight back to the original. Dito for obj.
  svref = vref;
  ssvref = vref; 
  svobj = vobj;
  // Initialise objects that will be used during the minimisation
  // robj is used for resampling vobj. Set size to that of ssvref
  robj = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(RefSz_x(),RefSz_y(),RefSz_z()));
  *robj = *ssvref;   // For header info
  *robj = 0.0;
  // datamask is a mask returned by general_transform with 1 when a voxel fall within the FOV and 0 otherwise
  datamask = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(RefSz_x(),RefSz_y(),RefSz_z()));
  copybasicproperties(*ssvref,*datamask);
  // df is a volume representation of the field
  df = boost::shared_ptr<NEWIMAGE::volume4D<float> >(new NEWIMAGE::volume4D<float>(RefSz_x(),RefSz_y(),RefSz_z(),df_field.size()));
  for (unsigned int i=0; i<df_field.size(); i++) {
    (*df)[i] = *ssvref;  // For header info
    (*df)[i] = 0.0;
    df_field[i]->AsVolume((*df)[i]);
  }
  // Various optional fields that have not yet been defined
  ref_fwhm = 0.0;                               // No smoothing yet
  obj_fwhm = 0.0;                               // No smoothing yet
  subsamp = std::vector<unsigned int>(3,1);     // No sub-sampling yet
  latest_ssd = 0.0;                             // No history
  lambda = 0.0;                                 // lambda not set
  mpl_lambda = 0.0;                             // Lambda for landmarks not set.
  use_ref_derivs = false;                       // Use "exact" derivatives as default
  hess_prec = BFMatrixDoublePrecision;          // Represent Hessian in double precision
  verbose = false;                              // Don't volunteer information
  debug = 0;                                    // Don't write debug info unless explicitly told to
  level = iter = attempt = 0;                   // Initilise debug info state variables
  regmod = BendingEnergy;                       // Regularise with bending energy as default
  ssd_lambda = false;                           // Don't weight lambda by ssd as default
  robj_updated = false;                         // Need to recalculate Robj
  robjmask_updated = false;                     // Need to recalculate Robj mask
  totmask_updated = false;                      // Need to recalculate Total mask
  robj_deriv_updated = false;                   // Need to recalculate derivatives of Robj
}

void fnirt_CF::SaveLocalIntensityMapping(std::string fname) const
{
  if (Verbose()) cout << "Saving Local Intensity Mapping to " << fname << endl;
  imap->SaveLocalMapping(fname);
}

void fnirt_CF::SaveGlobalIntensityMapping(std::string fname) const
{
  if (Verbose()) cout << "Saving Global Intensity Mapping to " << fname << endl;
  imap->SaveGlobalMapping(fname);
}

void fnirt_CF::SaveDefFields(std::string fname) const
{
  if (Verbose()) cout << "Saving Displacement fields to " << fname << endl;

  NEWIMAGE::FnirtFileWriter  write_it(fname,Ref(),DefFieldAsVol(0),DefFieldAsVol(1),DefFieldAsVol(2),AffineMat()); 
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
// This saves out a field of Jacobian-determinant values. It saves _only_ the Jacobian of 
// the non-linear part of the field. If you want to include also the affine part you need
// to calculate the determinant of the flirt matrix and add that to the values in the field.
//
////////////////////////////////////////////////////////////////////////////////////////////////

void fnirt_CF::SaveJacobian(std::string fname) const
{
  NEWIMAGE::volume<float>   jac(RefSz_x(),RefSz_y(),RefSz_z());
  jac = Ref();        // To get header-info
  jac = 0.0;          // Blank it
  deffield2jacobian(DefField(0),DefField(1),DefField(2),jac);

  if (Verbose()) cout << "Saving Jacobian map to " << fname << endl;
  save_volume(jac,fname);
}

void fnirt_CF::SaveDefCoefs(std::string fname) const
{
  if (Verbose()) cout << "Saving Coefficient fields to " << fname << endl;

  NEWIMAGE::FnirtFileWriter  write_it(fname,DefField(0),DefField(1),DefField(2),AffineMat());
}

void fnirt_CF::SaveRobj(std::string fname) const
{
  if (Verbose()) cout << "Saving Warped --in file to " << fname << endl;
  save_volume(Robj(),fname);
}

void fnirt_CF::SaveScaledRef(std::string fname) const
{
  NEWIMAGE::volume<float>        sref = Ref();  // To get all header info
  sref = 0.0;                                   // Blank it
  ScaledRef(sref);                              // Fill with scaled ref values
  if (Verbose()) cout << "Saving scaled --ref file to " << fname << endl;
  save_volume(sref,fname);
}
	     
void fnirt_CF::SaveRef(std::string fname) const
{
  if (Verbose()) cout << "Saving --ref file to " << fname << endl;
  save_volume(Ref(),fname);
}

void fnirt_CF::SaveMask(std::string fname) const
{
  if (Verbose()) cout << "Saving mask in --ref space to " << fname << endl;
  save_volume(Mask(),fname);
}

void fnirt_CF::SaveDataMask(std::string fname) const
{
  if (Verbose()) cout << "Saving data-mask in --ref space to " << fname << endl;
  save_volume(*datamask,fname);
}

void fnirt_CF::SaveRobjMask(std::string fname) const
{
  if (robjmask) {
    if (Verbose()) cout << "Saving warped --inmask in --ref space to " << fname << endl;
    save_volume(*robjmask,fname);
  }
}

void fnirt_CF::SaveObjMask(std::string fname) const
{
  if (objmask) {
    if (Verbose()) cout << "Saving --inmask to " << fname << endl;
    save_volume(*objmask,fname);
  }
}

std::pair<double,double> fnirt_CF::JacobianRange() const
{
  NEWIMAGE::volume<float>   jac(RefSz_x(),RefSz_y(),RefSz_z());
  jac = Ref();        // To get header-info
  jac = 0.0;          // Blank it
  deffield2jacobian(DefField(0),DefField(1),DefField(2),jac);
  std::pair<double,double>  rng;
  rng.first = jac.min();
  rng.second = jac.max();

  return(rng);
}

void fnirt_CF::ForceJacobianRange(double minj, double maxj) const
{
  NEWIMAGE::volume4D<float>   defvol(RefSz_x(),RefSz_y(),RefSz_z(),3);
  defvol.setdims(RefVxs_x(),RefVxs_y(),RefVxs_z(),1);
  for (unsigned int i=0; i<3; i++) defvol[i] = (*df)[i];
  // Use MJ's routines to constrain range of Jacobians.
  // N.B. that the range is not "guaranteed".
  convertwarp_rel2abs(defvol);
  constrain_topology(defvol,minj,maxj);
  convertwarp_abs2rel(defvol);
  // Set the constrained field as our "new" field
  SetDefField(defvol);
}

void fnirt_CF::SubsampleRef(const std::vector<unsigned int>& ss)
{
  if (ss != subsamp) { // If anything has really changed

    // Resample reference image and reference mask
    subsample_ref(ss2ms(ss),ss2vxs(ss));
    subsample_refmask(ss2ms(ss),ss2vxs(ss));

    // Inform user if requested
    if (Verbose()) {
      std::vector<unsigned int> new_ms = ss2ms(ss);
      std::vector<double> new_vxs = ss2vxs(ss);
      cout << "New Matrix Size: " << new_ms[0] << " " << new_ms[1] << " " << new_ms[2] << endl;
      cout << "New Voxel Size: " << new_vxs[0] << " " << new_vxs[1] << " " << new_vxs[2] << endl;
    }

    // Resample field
    for (unsigned int i=0; i<df_field.size(); i++) {
      boost::shared_ptr<BASISFIELD::basisfield> tmp = df_field[i]->ZoomField(ss2ms(ss),ss2vxs(ss)); 
      df_field[i] = tmp;
    }

    // Inform intensity-mapper of new dimensions
    imap->NewSubSampling(ss2ms(ss),ss2vxs(ss),ss,subsamp);

    // Resize volume representation of field
    df = boost::shared_ptr<NEWIMAGE::volume4D<float> >(new NEWIMAGE::volume4D<float>(RefSz_x(),RefSz_y(),RefSz_z(),df_field.size()));
    for (unsigned int i=0; i<df_field.size(); i++) {
      copybasicproperties(*ssvref,(*df)[i]);
      df_field[i]->AsVolume((*df)[i]);
    }    

    // Resize robj, robjmask and datamask
    robj = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(RefSz_x(),RefSz_y(),RefSz_z()));
    robj->copyproperties(*ssvref);
    robj_updated = false;
    robjmask = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(RefSz_x(),RefSz_y(),RefSz_z()));
    copybasicproperties(*ssvref,*robjmask);
    robjmask_updated = false;
    datamask = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(RefSz_x(),RefSz_y(),RefSz_z()));
    copybasicproperties(*ssvref,*datamask);

    // Free up memory from old derivatives
    if (robj_deriv) {
      robj_deriv = boost::shared_ptr<NEWIMAGE::volume4D<float> >();
      robj_deriv_updated = false;
    }    

    // Resize totmask
    if (totmask) {
      totmask = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(RefSz_x(),RefSz_y(),RefSz_z()));
      copybasicproperties(*ssvref,*totmask);
      totmask_updated = false;
    }

    subsamp = ss;
  }
}

void fnirt_CF::SmoothRef(double fwhm)
{
  if (fwhm != ref_fwhm) {
    ref_fwhm = fwhm;
    if (!fwhm) {
      svref = vref;
    }
    else {
      svref = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(smooth(*vref,fwhm/sqrt(8.0*log(2.0)))));
    }
    if (subsamp[0]!=1 || subsamp[1]!=1 || subsamp[2]!=1) {
      // Resample reference image and reference mask
      subsample_ref(ss2ms(subsamp),ss2vxs(subsamp));
      subsample_refmask(ss2ms(subsamp),ss2vxs(subsamp));
    }
    else ssvref = svref;
  }
  if (Verbose()) cout << "New FWHM (mm) for --ref: " << fwhm << endl;
}

void fnirt_CF::IntensityScaleObj(double sfac)
{
  *vobj *= sfac;      // Scale

  if (obj_fwhm) {     // If there is a smoothed version, smooth scaled image
    svobj = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(smooth(*vobj,obj_fwhm/sqrt(8.0*log(2.0)))));
  }
  robj_updated = false;
}

void fnirt_CF::SmoothObj(double fwhm)
{
  if (fwhm != obj_fwhm) {
    obj_fwhm = fwhm;
    if (!fwhm) {
      svobj = vobj;
    }
    else {
      svobj = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(smooth(*vobj,fwhm/sqrt(8.0*log(2.0)))));
    }
    robj_updated = false;
  }    

  if (Verbose()) cout << "New FWHM (mm) for --in: " << fwhm << endl;
}

void fnirt_CF::SetRefMask(const NEWIMAGE::volume<char>& mask)
{
  const boost::shared_ptr<NEWIMAGE::volume<char> >  tmp = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(mask));
  SetRefMask(tmp);
}

void fnirt_CF::SetRefMask(const boost::shared_ptr<NEWIMAGE::volume<char> > mask)
{
  if ((mask->xsize() != int(OrigRefSz_x())) || (mask->ysize() != int(OrigRefSz_y())) || (mask->zsize() != int(OrigRefSz_z()))) {
    FnirtException("fnirt_CF::set_refmask: Mismatch between reference volume and reference mask");
  }
  refmask = mask;
  copybasicproperties(*vref,*refmask);
  if (subsamp[0]!=1 || subsamp[1]!=1 || subsamp[2]!=1) {
    subsample_refmask(ss2ms(subsamp),ss2vxs(subsamp));
  }
  else ssrefmask = refmask;    
  if (!totmask) {
    totmask = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(RefSz_x(),RefSz_y(),RefSz_z()));
    copybasicproperties(*ssvref,*totmask);
  }
  totmask_updated = false;
}

void fnirt_CF::SetObjMask(const NEWIMAGE::volume<char>& mask)
{
  const boost::shared_ptr<NEWIMAGE::volume<char> >  tmp = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(mask));
  SetObjMask(tmp);
}

void fnirt_CF::SetObjMask(const boost::shared_ptr<NEWIMAGE::volume<char> > mask)
{
  if ((mask->xsize() != int(ObjSz_x())) || (mask->ysize() != int(ObjSz_y())) || (mask->zsize() != int(ObjSz_z()))) {
    FnirtException("fnirt_CF::set_refmask: Mismatch between reference volume and reference mask");
  }
  objmask = mask;
  copybasicproperties(*vobj,*objmask);
  robjmask = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(RefSz_x(),RefSz_y(),RefSz_z()));
  copybasicproperties(*ssvref,*robjmask);
  robjmask_updated = false;
  if (!totmask) {
    totmask = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(RefSz_x(),RefSz_y(),RefSz_z()));
    copybasicproperties(*ssvref,*totmask);
  }
  totmask_updated = false;
}

NEWMAT::ColumnVector fnirt_CF::GetDefFieldParams(int  findx) const
{
  if (findx < -1 || findx > 2) {
    FnirtException("fnirt_CF::GetDefFieldParams: Field index out of range");
  }
  if (findx == -1) {
    NEWMAT::ColumnVector  ret;
    for (int i=0; i<3; i++) {
      ret &= *(DefField(i).GetCoef());
    }
    return(ret);
  }
  return(*(DefField(findx).GetCoef()));
}

NEWMAT::ColumnVector fnirt_CF::GetScaleParams() const
{
  NEWMAT::ColumnVector    scpar = imap->GetPar();
  return(scpar);
}

void fnirt_CF::update_robj() const
{
  if (!robj_updated) {
    general_transform(*svobj,aff,*df,*robj,*datamask);
    robj_updated = true;
    totmask_updated = false;
  }    
}

void fnirt_CF::update_robjmask() const
{
  if (objmask && !robjmask_updated) {
    general_transform(*objmask,aff,*df,*robjmask);
    robjmask_updated = true;
    totmask_updated = false;
  }
}

const NEWIMAGE::volume<float>& fnirt_CF::Robj() const
{
  update_robj();
  return(*robj);
}

const NEWIMAGE::volume4D<float>& fnirt_CF::RobjDeriv() const
{
  if (!robj_deriv) { // If first time derivative of object image is requested
    robj_deriv = boost::shared_ptr<NEWIMAGE::volume4D<float> >(new NEWIMAGE::volume4D<float>(RefSz_x(),RefSz_y(),RefSz_z(),3));
    for (int i=0; i<3; i++) (*robj_deriv)[i].copyproperties(*ssvref);
    robj_deriv_updated = false;
  }
  if (!robj_deriv_updated) { // If the derivatives doesn't reflect current df field
    // When we are doing this we might as well update robj 
    // and data mask as well.
    general_transform_3partial(*svobj,aff,*df,*robj,*robj_deriv,*datamask);
    (*robj_deriv)[0] /= ObjVxs_x();
    (*robj_deriv)[1] /= ObjVxs_y();
    (*robj_deriv)[2] /= ObjVxs_z();
    if (!robj_updated) {robj_updated=true; totmask_updated=false;}
    robj_deriv_updated = true;
  }
  return(*robj_deriv);
}

const NEWIMAGE::volume<float>& fnirt_CF::Ref() const
{
  return(*ssvref);
}

void fnirt_CF::ScaledRef(NEWIMAGE::volume<float>&  sref) const
{
  if (sref.xsize() != Ref().xsize() || sref.ysize() != Ref().ysize() || sref.zsize() != Ref().zsize() ||
      sref.xdim() != Ref().xdim() || sref.ydim() != Ref().ydim() || sref.zdim() != Ref().zdim() ) {
    throw FnirtException("fnirt_CF::ScaledRef: Incompatible size of sref");
  }
  imap->ScaleMe(Ref(),sref);
}

const NEWIMAGE::volume4D<float>& fnirt_CF::RefDeriv() const
{
  if (!ref_deriv) {
    ref_deriv = boost::shared_ptr<NEWIMAGE::volume4D<float> >(new NEWIMAGE::volume4D<float>(RefSz_x(),RefSz_y(),RefSz_z(),3));
    ref_deriv->copyproperties(*ssvref);
    NEWIMAGE::volume<float>    skrutt(RefSz_x(),RefSz_y(),RefSz_z());    
    skrutt.copyproperties(*ssvref);
    affine_transform_3partial(*ssvref,IdentityMatrix(4),skrutt,*ref_deriv);
    (*ref_deriv)[0] /= RefVxs_x();
    (*ref_deriv)[1] /= RefVxs_y();
    (*ref_deriv)[2] /= RefVxs_z();
  }
  return(*ref_deriv);
}

const NEWIMAGE::volume4D<float>& fnirt_CF::Deriv() const
{
  if (use_ref_derivs) {
    return(RefDeriv());
  }
  else {
    return(RobjDeriv());
  }
}

// Return a "total" mask from data mask, mask in reference space and mask in object space
const NEWIMAGE::volume<char>& fnirt_CF::Mask() const
{
  if (totmask) {
    if (!totmask_updated) {
      update_robj();      // Bring data mask up to date
      update_robjmask();  // Bring object mask up to date
      // Now combine the masks
      for (unsigned int k=0; k<RefSz_z(); k++) {
        for (unsigned int j=0; j<RefSz_y(); j++) {
          for (unsigned int i=0; i<RefSz_x(); i++) {
            (*totmask)(i,j,k) = (ssrefmask) ? (*datamask)(i,j,k) & (*ssrefmask)(i,j,k) : (*datamask)(i,j,k);
            if (objmask) {(*totmask)(i,j,k) = ((*robjmask)(i,j,k)>0.5) ? (*totmask)(i,j,k) : 0;}
          }
        }
      }
      totmask_updated = true;
    }
    return(*totmask);
  }
  // else 
  return(*datamask);
}

void fnirt_CF::SetDefField(const NEWIMAGE::volume4D<float>& vfield) const
{
  if (vfield.xsize() != int(RefSz_x()) || vfield.ysize() != int(RefSz_y()) || vfield.zsize() != int(RefSz_z())) {
    throw FnirtException("fnirt_CF::SetDefField: Size mismatch between current field field and parameter vector");
  }
  if (fabs(vfield.xdim()-RefVxs_x()) > 1e-6 || fabs(vfield.ydim()-RefVxs_y()) > 1e-6 || fabs(vfield.zdim()-RefVxs_z()) > 1e-6) {
    throw FnirtException("fnirt_CF::SetDefField: Voxel-size mismatch between current field field and parameter vector");
  }
  // Might want to introduce check for equality in future
  for (unsigned int i=0; i<df_field.size(); i++) {
    df_field[i]->Set(vfield[i]);
    df_field[i]->AsVolume((*df)[i]);   // N.B. may not be identical to vfield[i]
  }
  robj_updated = false;
  robj_deriv_updated = false;
  if (objmask) {
    totmask_updated = false;
    robjmask_updated = false;
  }
}

void fnirt_CF::SetDefFieldParams(const NEWMAT::ColumnVector& p) const
{
  if (p.Nrows() != 3*int(DefCoefSz())) {throw FnirtException("fnirt_CF::SetDefFieldParams: Mismatch between field and parameter vector");}
  bool sameasold = true;
  for (unsigned int i=0; i<df_field.size(); i++) {
    boost::shared_ptr<NEWMAT::ColumnVector>  oldp = df_field[i]->GetCoef();  // Shared ptr must be named variable
    if (*oldp != p.Rows(i*DefCoefSz()+1,(i+1)*DefCoefSz())) {sameasold = false;}
  }

  if (!sameasold) {
    for (int i=0; i<int(df_field.size()); i++) {
      df_field[i]->SetCoef(p.Rows(i*DefCoefSz()+1,(i+1)*DefCoefSz()));
      df_field[i]->AsVolume((*df)[i]);
    }
    robj_updated = false;
    robj_deriv_updated = false;
    if (objmask) {
      totmask_updated = false;
      robjmask_updated = false;
    }
  }  
}

void fnirt_CF::SetScaleParams(const NEWMAT::ColumnVector& p)
const
{
  imap->SetPar(p);
}

void fnirt_CF::subsample_ref(const std::vector<unsigned int>  nms,
                             const std::vector<double>        nvxs)
{
  ssvref = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(int(nms[0]),int(nms[1]),int(nms[2])));
  ssvref->copyproperties(*svref);
  ssvref->setdims(nvxs[0],nvxs[1],nvxs[2]);

  // Set sform and qform appropriately. What we want are sform and qform
  // matrices that have the new voxel-sizes, and that maps the [0 0 0] point
  // To exactly the same point (as before) in world space.
  NEWMAT::Matrix M(4,4);
  M  = 0.0;
  M(1,1) = nvxs[0]/svref->xdim(); M(2,2) = nvxs[1]/svref->ydim(); M(3,3) = nvxs[2]/svref->zdim(); M(4,4) = 1.0;
  if (svref->sform_code() != NIFTI_XFORM_UNKNOWN) ssvref->set_sform(svref->sform_code(),svref->sform_mat() * M);
  if (svref->qform_code() != NIFTI_XFORM_UNKNOWN) ssvref->set_qform(svref->qform_code(),svref->qform_mat() * M);

  // And now resample. Note that we don't do any integration, but rather assume
  // that the smoothing of the original volume above has effectively done that.
  NEWIMAGE::volume4D<float>  fakevol;
  std::vector<int>           fakevec;
  raw_general_transform(*svref,IdentityMatrix(4),fakevol,fakevec,fakevec,*ssvref,fakevol);
}

void fnirt_CF::subsample_refmask(const std::vector<unsigned int>  nms,
                                 const std::vector<double>        nvxs)
{
  if (refmask) {
    ssrefmask = boost::shared_ptr<NEWIMAGE::volume<char> >(new NEWIMAGE::volume<char>(int(nms[0]),int(nms[1]),int(nms[2])));
    ssrefmask->copyproperties(*refmask);
    ssrefmask->setdims(nvxs[0],nvxs[1],nvxs[2]);

    // Set sform and qform appropriately. What we want are sform and qform
    // matrices that have the new voxel-sizes, and that maps the [0 0 0] point
    // To exactly the same point (as before) in world space.
    NEWMAT::Matrix M(4,4);
    M  = 0.0;
    M(1,1) = nvxs[0]/refmask->xdim(); M(2,2) = nvxs[1]/refmask->ydim(); M(3,3) = nvxs[2]/refmask->zdim(); M(4,4) = 1.0;
    if (refmask->sform_code() != NIFTI_XFORM_UNKNOWN) ssrefmask->set_sform(refmask->sform_code(),refmask->sform_mat() * M);
    if (refmask->qform_code() != NIFTI_XFORM_UNKNOWN) ssrefmask->set_qform(refmask->qform_code(),refmask->qform_mat() * M);

    // I will (rather lazily) use a temporary float-copy of the mask for the 
    // downsampling. This will then be thresholded at 0.99 when translating back to char.
    NEWIMAGE::volume<float> finmask(vref->xsize(),vref->ysize(),vref->zsize());
    copybasicproperties(*refmask,finmask);
    finmask.insert_vec(refmask->vec());
    NEWIMAGE::volume<float> foutmask(ssrefmask->xsize(),ssrefmask->ysize(),ssrefmask->zsize());
    copybasicproperties(*ssrefmask,foutmask);
    
    // And now resample. Note that we don't do any integration, but rather assume
    // that the smoothing of the original volume above has effectively done that.
    // Note also that after this call ssrefmask will be one for all voxels falling
    // inside the original volume and zero otherwise.
    NEWIMAGE::volume4D<float>  fakevol;
    std::vector<int>           fakevec;
    raw_general_transform(finmask,IdentityMatrix(4),fakevol,fakevec,fakevec,foutmask,fakevol,*ssrefmask);

    // Now combine the resampled mask with the "inside volume" mask;
    for (int k=0; k<ssrefmask->zsize(); k++) {
      for (int j=0; j<ssrefmask->ysize(); j++) {
	for (int i=0; i<ssrefmask->xsize(); i++) {
          (*ssrefmask)(i,j,k) = (foutmask(i,j,k) > 0.99) ? (*ssrefmask)(i,j,k) : 0;
	}
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
// Member functions for class SSD_fnirt_CF
//
///////////////////////////////////////////////////////////////////////////////////////////////

NEWMAT::ReturnMatrix SSD_fnirt_CF::Par() const
{
  NEWMAT::ColumnVector ovec(NPar());

  ovec.Rows(1,3*DefCoefSz()) = GetDefFieldParams();
  ovec.Rows(3*DefCoefSz()+1,3*DefCoefSz()+ScaleCoefSz()) = IMap()->GetPar();

  ovec.Release();
  return(ovec);
}

double SSD_fnirt_CF::cf(const NEWMAT::ColumnVector& p) const
{
  if (p.Nrows() != NPar()) {
    throw FnirtException("SSD_fnirt_CF::cf: Mismatch between field and parameter vector");
  }
  SetDefFieldParams(p.Rows(1,3*DefCoefSz()));
  if (!IMap()->Fixed()) {
    SetScaleParams(p.Rows(3*DefCoefSz()+1,3*DefCoefSz()+ScaleCoefSz()));
  }
  const NEWIMAGE::volume<float>& ref = Ref();
  const NEWIMAGE::volume<float>& obj = Robj();
  const NEWIMAGE::volume<char>&  mask = Mask();
  
  NEWIMAGE::volume<float>        sref(ref.xsize(),ref.ysize(),ref.zsize());
  copybasicproperties(ref,sref);
  ScaledRef(sref);                              // Returns intensity-mapped ref in sref.

  double ssd = 0.0;
  int n = 0;
  for (unsigned int k=0; k<RefSz_z(); k++) {
    for (unsigned int j=0; j<RefSz_y(); j++) {
      for (unsigned int i=0; i<RefSz_x(); i++) {
        if (mask(i,j,k)) {
          n++; 
          ssd += SQR(obj(i,j,k) - sref(i,j,k));
        }
      }
    }
  }
  ssd /= double(n);                             // Mean SSD
  SetLatestSSD(ssd);

  if (Debug()) {
    SetAttempt(Attempt()+1);
    bool old_verbose = Verbose();
    LocalSetVerbose(false);
    SaveRobj(string("FnirtDebugRobj_")+DebugString());
    SaveScaledRef(string("FnirtDebugScaledRef_")+DebugString());
    SaveJacobian(string("FnirtDebugJacobian_")+DebugString());
    if (Debug() > 1) {    
      SaveMask(string("FnirtDebugMask_")+DebugString());
      SaveDataMask(string("FnirtDebugDataMask_")+DebugString());
      SaveRobjMask(string("FnirtDebugRobjMask_")+DebugString());
    }
    LocalSetVerbose(old_verbose);
  }

  double cost = ssd;
  if (Verbose()) cout << "SSD = " << ssd;

  double disp_reg = 0.0;
  for (int i=0; i<3; i++) {
    if (RegularisationModel() == MembraneEnergy) {
      disp_reg += Lambda()*DefField(i).MemEnergy()/double(n);  // "Average" membrane energy
    }
    else if (RegularisationModel() == BendingEnergy) {
      disp_reg += Lambda()*DefField(i).BendEnergy()/double(n);  // "Average" bending energy
    }
  }
  if (Verbose()) cout << "\tDisp reg = " << disp_reg;
  cost += disp_reg;

  double int_reg = IMap()->CFContribution() / double(n);           // Add regularisation of intensity mapping
  if (Verbose()) cout << "\tInt reg = " << int_reg;
  cost += int_reg;

  if (MPL()) {                                                    // Add (optional) squared landmark distance
    double mp_ssd = 0.0;
    for (unsigned int i=0; i<3; i++) mp_ssd += MatchingPointsLambda() * MPL()->SSD(i,DefField(i));
    if (Verbose()) cout << "\tLM-dist = " << mp_ssd;
    cost += mp_ssd;
  }
  if (Verbose()) cout << "\tTotal Cost = " << cost << endl;

  return(cost);
}

// Gradient of cost function

NEWMAT::ReturnMatrix SSD_fnirt_CF::grad(const NEWMAT::ColumnVector& p) const
{
  if (p.Nrows() != NPar()) {
    throw FnirtException("SSD_fnirt_CF::cf: Mismatch between field and parameter vector");
  }
  SetDefFieldParams(p.Rows(1,3*DefCoefSz()));
  if (!IMap()->Fixed()) {
    SetScaleParams(p.Rows(3*DefCoefSz()+1,3*DefCoefSz()+ScaleCoefSz()));
  }
  NEWMAT::ColumnVector  gradient(NPar());

  const NEWIMAGE::volume<float>& ref = Ref();
  const NEWIMAGE::volume<float>& obj = Robj();
  const NEWIMAGE::volume<char>&  mask = Mask();
  unsigned int n = static_cast<unsigned int>(mask.sum());
  NEWIMAGE::volume<float>        sref(ref.xsize(),ref.ysize(),ref.zsize());
  copybasicproperties(ref,sref);
  ScaledRef(sref);               // Returns intensity-mapped ref in sref

  const NEWIMAGE::volume4D<float>& partial = Deriv();    // Calls either RobjDeriv() or RefDeriv()
  for (int v=0; v<3; v++) {
    gradient.Rows(v*DefCoefSz()+1,(v+1)*DefCoefSz()) = 2.0 * DefField(v).Jte(partial[v],obj-sref,&mask) / double(n);         // Gradient of mean-SSD

    // Add regularisation component
    if (RegularisationModel() == MembraneEnergy) {
      gradient.Rows(v*DefCoefSz()+1,(v+1)*DefCoefSz()) += Lambda()*DefField(v).MemEnergyGrad() / double(n); // Gradient of "average" membrane energy
    }
    else if (RegularisationModel() == BendingEnergy) {
      gradient.Rows(v*DefCoefSz()+1,(v+1)*DefCoefSz()) += Lambda()*DefField(v).BendEnergyGrad() / double(n); // Gradient of "average" bending energy
    }

    // Add (optional) point-list/landmarks component
    if (MPL()) {
      // cout << "Adding bit from point lists. MPLambda = " << MatchingPointsLambda() << endl;
      gradient.Rows(v*DefCoefSz()+1,(v+1)*DefCoefSz()) += MatchingPointsLambda() * MPL()->SSD_Gradient(static_cast<unsigned int>(v),DefField(v));
    }
  }

  // Bottom-part of gradient is derivative w.r.t. scaling parameters

  NEWMAT::ColumnVector scgrad = IMap()->Gradient(ref,obj-sref,&mask);
  gradient.Rows(3*DefCoefSz()+1,3*DefCoefSz()+ScaleCoefSz()) = (2.0/double(n)) * scgrad;

  if (Debug()) {
    SetIter(Iter()+1);
    if (Debug() > 1) {
      MISCMATHS::write_ascii_matrix(string("FnirtDebugGradient_")+DebugString()+string(".txt"),gradient);
    }
  }
  
  gradient.Release();
  return(gradient);
}

// Gauss-Newton approximation to hessian of cost function
// The hessian is built in pieces as
// [(df/dx)^T(df/dx) (df/dx)^T(df/dy) (df/dx)^T(df/dz)]
// [(df/dy)^T(df/dx) (df/dy)^T(df/dy) (df/dy)^T(df/dz)]
// [(df/dz)^T(df/dx) (df/dz)^T(df/dy) (df/dz)^T(df/dz)]

boost::shared_ptr<MISCMATHS::BFMatrix> SSD_fnirt_CF::hess(const NEWMAT::ColumnVector&             p,
                                                          boost::shared_ptr<MISCMATHS::BFMatrix>  iptr) const
{
  //
  // See if there is an old Hessian, and if we can reuse it
  //
  if (last_hess && UsingRefDeriv()) {               // If we are using derivatives from reference scan
    if (int(last_hess->Ncols()) == p.Nrows()) {     // If there has been no zooming
      return(last_hess);
    }
    else { 
      last_hess->Clear();
      last_hess = boost::shared_ptr<MISCMATHS::BFMatrix>();
    }
  }
  // Reclaim memory from last Hessian
  if (iptr) {
    iptr->Clear();
  }

  if (p.Nrows() != NPar()) {
    throw FnirtException("SSD_fnirt_CF::hess: Mismatch between field and parameter vector");
  }
  //
  // This means we have to calculate it. Phew!
  //
  SetDefFieldParams(p.Rows(1,3*DefCoefSz()));
  if (!IMap()->Fixed()) {
    SetScaleParams(p.Rows(3*DefCoefSz()+1,3*DefCoefSz()+ScaleCoefSz()));
  }
  //
  // Get masked ref, and derivatives in Column vectors.
  //
  const NEWIMAGE::volume<float>& ref = Ref();
  const NEWIMAGE::volume4D<float>& partial = Deriv();    // Calls RobjDeriv or RefDeriv
  const NEWIMAGE::volume<char>& mask = Mask();
  unsigned int n = static_cast<unsigned int>(mask.sum());

  //
  //Start building H
  //
  // Diagonal blocks
  boost::shared_ptr<MISCMATHS::BFMatrix> dxTdx = DefField(0).JtJ(partial[0],&mask,HessianPrecision());
  boost::shared_ptr<MISCMATHS::BFMatrix> dyTdy = DefField(0).JtJ(partial[1],&mask,HessianPrecision());
  boost::shared_ptr<MISCMATHS::BFMatrix> dzTdz = DefField(0).JtJ(partial[2],&mask,HessianPrecision());
  dxTdx->MulMeByScalar(2.0/double(n));
  dyTdy->MulMeByScalar(2.0/double(n));
  dzTdz->MulMeByScalar(2.0/double(n));

  // Off-diagonal blocks
  boost::shared_ptr<MISCMATHS::BFMatrix> dxTdy = DefField(0).JtJ(partial[1],partial[0],&mask,HessianPrecision());
  boost::shared_ptr<MISCMATHS::BFMatrix> dxTdz = DefField(0).JtJ(partial[2],partial[0],&mask,HessianPrecision());
  boost::shared_ptr<MISCMATHS::BFMatrix> dyTdz = DefField(0).JtJ(partial[2],partial[1],&mask,HessianPrecision());
  dxTdy->MulMeByScalar(2.0/double(n));
  dxTdz->MulMeByScalar(2.0/double(n));
  dyTdz->MulMeByScalar(2.0/double(n));

  // Add regularisation to diagonal blocks
  if (!reg || reg->Ncols() != dxTdx->Ncols()) {
    if (reg) reg->Clear();
    if (RegularisationModel() == MembraneEnergy) reg = DefField(0).MemEnergyHess(HessianPrecision());
    else if (RegularisationModel() == BendingEnergy) reg = DefField(0).BendEnergyHess(HessianPrecision());
    if (Debug() > 1) reg->Print(string("FnirtDebugRegHess_")+DebugString()+string(".txt"));
  }
  dxTdx->AddToMe(*reg,Lambda()/double(n));
  dyTdy->AddToMe(*reg,Lambda()/double(n));
  dzTdz->AddToMe(*reg,Lambda()/double(n));

  // Add (optional) contribution from point-list/landmarks.
  if (MPL()) {
    boost::shared_ptr<MISCMATHS::BFMatrix> mplH = MPL()->SSD_Hessian(0,DefField(0),HessianPrecision());
    dxTdx->AddToMe(*mplH,MatchingPointsLambda());
    if (Debug() > 1) mplH->Print(string("FnirtDebugMPL_x_Hess_")+DebugString()+string(".txt"));
    mplH = MPL()->SSD_Hessian(1,DefField(1),HessianPrecision());
    dyTdy->AddToMe(*mplH,MatchingPointsLambda());
    if (Debug() > 1) mplH->Print(string("FnirtDebugMPL_y_Hess_")+DebugString()+string(".txt"));
    mplH = MPL()->SSD_Hessian(2,DefField(2),HessianPrecision());
    dzTdz->AddToMe(*mplH,MatchingPointsLambda());
    if (Debug() > 1) mplH->Print(string("FnirtDebugMPL_z_Hess_")+DebugString()+string(".txt"));
  }

  // Put all the bits together, and lets put them all in dxTdx.
  dxTdx->HorConcat2MyRight(*dxTdy);
  dxTdx->HorConcat2MyRight(*dxTdz);              // dxTdx now contains upper three submatrices
  dxTdy->HorConcat2MyRight(*dyTdy);
  dyTdy->Clear();
  dxTdy->HorConcat2MyRight(*dyTdz);              // dxTdy now contains middle three submatrices
  dxTdz->HorConcat2MyRight(*dyTdz);
  dyTdz->Clear();
  dxTdz->HorConcat2MyRight(*dzTdz);              // DxTdz now contains lower three submatrices
  dzTdz->Clear();
  dxTdx->VertConcatBelowMe(*dxTdy);
  dxTdy->Clear();
  dxTdx->VertConcatBelowMe(*dxTdz);              // dxTdx now contains all 3x3 submatrices
  dxTdz->Clear();

  // Calculate part pertaining to intensity mapping.
  boost::shared_ptr<MISCMATHS::BFMatrix> imapH = IMap()->Hessian(ref,&mask,HessianPrecision());
  imapH->MulMeByScalar(2.0/double(n));

  // Calculate Cross-term between spatial and intensity mapping.
  boost::shared_ptr<MISCMATHS::BFMatrix> xCross = IMap()->CrossHessian(DefField(0),partial[0],ref,&mask,HessianPrecision());  
  xCross->MulMeByScalar(2.0/double(n));
  boost::shared_ptr<MISCMATHS::BFMatrix> yCross = IMap()->CrossHessian(DefField(0),partial[1],ref,&mask,HessianPrecision());  
  yCross->MulMeByScalar(2.0/double(n));
  boost::shared_ptr<MISCMATHS::BFMatrix> zCross = IMap()->CrossHessian(DefField(0),partial[2],ref,&mask,HessianPrecision());
  zCross->MulMeByScalar(2.0/double(n));

  // Collect them in xCross
  xCross->VertConcatBelowMe(*yCross);
  yCross->Clear();
  xCross->VertConcatBelowMe(*zCross);
  zCross->Clear();

  // And put it all together
  dxTdx->VertConcatBelowMe(*(xCross->Transpose()));
  xCross->VertConcatBelowMe(*imapH);
  dxTdx->HorConcat2MyRight(*xCross);
  xCross->Clear();
    
  if (UsingRefDeriv()) last_hess = dxTdx;

  if (Debug() > 1) dxTdx->Print(string("FnirtDebugHessian_")+DebugString()+string(".txt"));

  return(dxTdx);
}


