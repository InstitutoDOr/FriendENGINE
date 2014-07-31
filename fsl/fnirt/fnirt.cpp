// FNIRT - FMRIB's Non-linear Image Registration Tool
//
// fnirt.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
//

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/nonlin.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "matching_points.h"
#include "fnirt_costfunctions.h"
#include "intensity_mappers.h"
#include "fnirtfns.h"
#include "parser.h"

using namespace std;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace BASISFIELD;
using namespace FNIRT;


extern "C" __declspec(dllexport) int _stdcall fnirt(char *CmdLn)
{
  int argc;
  char **argv;
  
  parser(CmdLn, argc, argv);
  // Parse command line and return object that we can 
  // later ask what the user wants us to do.
  boost::shared_ptr<fnirt_clp> clp;
  try {
    clp = parse_fnirt_command_line(argc,argv);
  }
  catch (const std::exception& error) {
    cerr << "Error occured when parsing the command line" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    return(EXIT_FAILURE);
  }

  // Read --ref and --in volumes
  boost::shared_ptr<volume<float> >     ref(new volume<float>);
  boost::shared_ptr<volume<float> >     obj(new volume<float>);
  try {
    read_volume(*ref,clp->Ref());
    read_volume(*obj,clp->Obj());
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when reading --ref or --obj file" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    return(EXIT_FAILURE);
  }

  // Check for registration to self (sometimes done in TBSS)
  if (trying_to_register_to_self(clp->Ref(),*ref,clp->Obj(),*obj,clp->Affine())) {
    try {
      write_self_results(*clp,*ref);
    } 
    catch(const std::exception& error) {
      cerr << "Error occurred when writing results of registration to self" << endl;
      cerr << "Exception thrown with message: " << error.what() << endl; 
    }
    if (clp->Verbose()) cerr << "Input to fnirt was two identical images. Result is unity transform" << endl;
    return(EXIT_SUCCESS);  // Sort of
  }

  // Prepare for fnirting
  boost::shared_ptr<SSD_fnirt_CF>       cf;
  boost::shared_ptr<NonlinParam>        nlpar;
  double                                objmean = 0.0;
  try {
    // Create refmask as combination of explicit and implicit masks.
    boost::shared_ptr<volume<char> > refmask;
    if (clp->UseRefMask(1)) refmask = make_mask(clp->RefMask(),InclusiveMask,*ref,clp->UseImplicitRefMask(),clp->ImplicitRefValue());
    else refmask = make_mask(clp->RefMask(),IgnoreMask,*ref,clp->UseImplicitRefMask(),clp->ImplicitRefValue());

    // Normalise ref global mean to 100
    double refmean = spmlike_mean(*ref);
    (*ref) *= (100.0/refmean);     

    // Create objmask as a combination of explicit and implicit masks.
    boost::shared_ptr<volume<char> > objmask;
    if (clp->UseObjMask(1)) objmask = make_mask(clp->ObjMask(),InclusiveMask,*obj,clp->UseImplicitObjMask(),clp->ImplicitObjValue());
    else objmask = make_mask(clp->ObjMask(),IgnoreMask,*obj,clp->UseImplicitObjMask(),clp->ImplicitObjValue());

    // Normalise obj to global mean 100
    objmean = spmlike_mean(*obj);
    (*obj) *= (100.0/objmean);          

    // Initialise the field that describes the warps.
    std::vector<boost::shared_ptr<basisfield> >    field = init_warpfield(*clp);

    // Set up model for intensity mapping
    boost::shared_ptr<IntensityMapper>   mymapper = init_intensity_mapper(*clp);

    // Create cost-function object
    cf = boost::shared_ptr<SSD_fnirt_CF>(new SSD_fnirt_CF(*ref,*obj,clp->Affine(),field,mymapper));
    cf->SetVerbose(clp->Verbose());
    cf->SetRegularisationModel(clp->RegularisationModel());
    cf->SetLambda(clp->Lambda(1));
    cf->SetLevel(1);
    cf->SetIntensityMappingFixed(!clp->EstimateIntensity(1));
    cf->SetHessianPrecision(clp->HessianPrecision());
    if (clp->WeightLambdaBySSD()) cf->WeightLambdaBySSD();
    if (clp->UseRefDeriv()) cf->UseRefDerivs();
    if (clp->Debug()) cf->SetDebug(clp->Debug());

    // Smooth and sub-sample
    cf->SmoothRef(clp->RefFWHM(1));
    cf->SmoothObj(clp->ObjFWHM(1));
    cf->SubsampleRef(clp->SubSampling(1));

    // Set masks (if any)
    if (refmask) cf->SetRefMask(*refmask);
    if (objmask) cf->SetObjMask(*objmask);

    // Check if we should include point-lists in the warping
    if (clp->RefPointList().length()) {
      PointList         rpl(clp->RefPointList(),clp->Ref());
      rpl.SetAffine(clp->Affine().i());
      PointList         opl(clp->ObjPointList(),clp->Obj());
      MatchingPoints    mpl(rpl,opl);
      cf->SetMatchingPoints(mpl);
      cf->SetMatchingPointsLambda(clp->MatchingPointsLambda(1));
    }

    // Initialise nonlin parameter object
    NEWMAT::ColumnVector             spar = cf->Par();   // Start guesses for parameters
    nlpar = boost::shared_ptr<NonlinParam>(new NonlinParam(cf->NPar(),clp->MinimisationMethod(),spar));
    set_nlpars(*nlpar);
    nlpar->SetMaxIter(clp->MaxIter(1));
  }
  catch (const std::exception& error) {
    cerr << "Error occurred when preparing to fnirt" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    return(EXIT_FAILURE);
  }

  // And now find parameters for the first level of subsampling
  try {
    // MISCMATHS::NonlinOut status = nonlin(*nlpar,*cf);
    nonlin(*nlpar,*cf);
    if (!constrain_warpfield(*cf,*clp,5)) { // N.B. arbitrary number 5
      std::pair<double,double> range = cf->JacobianRange();
      cout << "Warning, Jacobian not within prescribed range. Prescription is " << clp->JacLowerBound() << " -- " << clp->JacUpperBound();
      cout << " and obtained range is " << range.first << " -- " << range.second << endl;
    }
  }
  catch (const std::exception& error) {
    cerr << "Error occured during estimation at first level of subsampling" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    return(EXIT_FAILURE);
  }

  // Loop over remaining levels of sub-sampling
  try {
    for (int ssl=2; ssl<=clp->NoLevels(); ssl++) {
      if (clp->Verbose()) cout << "***Going to next resolution level***" << endl;
      cf->SetLevel(ssl);

      // New mask use
      if (clp->UseRefMask(ssl) != clp->UseRefMask(ssl-1)) {
        boost::shared_ptr<volume<char> > refmask;
        if (clp->UseRefMask(ssl)) refmask = make_mask(clp->RefMask(),InclusiveMask,*ref,clp->UseImplicitRefMask(),clp->ImplicitRefValue());
        else refmask = make_mask(clp->RefMask(),IgnoreMask,*ref,clp->UseImplicitRefMask(),clp->ImplicitRefValue());
        if (refmask) cf->SetRefMask(*refmask);      
      }
      if (clp->UseObjMask(ssl) != clp->UseObjMask(ssl-1)) {
        boost::shared_ptr<volume<char> > objmask;
        if (clp->UseObjMask(1)) objmask = make_mask(clp->ObjMask(),ExclusiveMask,*obj,clp->UseImplicitObjMask(),clp->ImplicitObjValue());
        else objmask = make_mask(clp->ObjMask(),IgnoreMask,*obj,clp->UseImplicitObjMask(),clp->ImplicitObjValue());
        if (objmask) cf->SetObjMask(*objmask);
      }

      // New smoothness
      if (clp->RefFWHM(ssl) != clp->RefFWHM(ssl-1)) cf->SmoothRef(clp->RefFWHM(ssl));
      if (clp->ObjFWHM(ssl) != clp->ObjFWHM(ssl-1)) cf->SmoothObj(clp->ObjFWHM(ssl));

      // New subsampling
      cout << "Setting subsampling" << endl;
      if (clp->SubSampling(ssl) != clp->SubSampling(ssl-1)) cf->SubsampleRef(clp->SubSampling(ssl));
      cout << "Setting reg mode" << endl;
      cf->SetRegularisationModel(clp->RegularisationModel());
      cout << "Setting lambda" << endl;
      cf->SetLambda(clp->Lambda(ssl));  // New (possibly) lambda
      cf->SetMatchingPointsLambda(clp->MatchingPointsLambda(ssl));
      cf->SetIntensityMappingFixed(!clp->EstimateIntensity(ssl));

      // Splash out on brand new nonlin object
      nlpar = boost::shared_ptr<NonlinParam>(new NonlinParam(cf->NPar(),clp->MinimisationMethod(),cf->Par()));
      set_nlpars(*nlpar);
      nlpar->SetMaxIter(clp->MaxIter(ssl));

      // MISCMATHS::NonlinOut status = nonlin(*nlpar,*cf);
      nonlin(*nlpar,*cf);
      if (!constrain_warpfield(*cf,*clp,5)) { // N.B. arbitrary number 5
        std::pair<double,double> range = cf->JacobianRange();
        cout << "Warning, Jacobian not within prescribed range. Prescription is " << clp->JacLowerBound() << " -- " << clp->JacUpperBound();
        cout << " and obtained range is " << range.first << " -- " << range.second << endl;
      }
    }

    // If we stopped short of full resolution, 
    // up-sample field to full resolution.
    if (clp->SubSampling(clp->NoLevels()) > 1) cf->SubsampleRef(1);
  }
  catch (const std::exception& error) {
    cerr << "Error occured during estimation at subsampling level > 1" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    return(EXIT_FAILURE);
  }

  // Save everything we have been asked to save
  try {
    cf->SaveDefCoefs(clp->CoefFname());                                      // Coefficients
    if (clp->FieldFname().length()) cf->SaveDefFields(clp->FieldFname());    // Field
    if (clp->JacFname().length()) cf->SaveJacobian(clp->JacFname());         // Jacobian
    if (clp->RefOutFname().length()) cf->SaveScaledRef(clp->RefOutFname());  // Intensity modulated ref scan
    // Intensity-mapping
    if (clp->IntensityMappingFname().length()) cf->SaveIntensityMapping(clp->IntensityMappingFname());
    // Warped object image
    if (clp->ObjOutFname().length()) {
      if (clp->ObjFWHM(clp->NoLevels())!=0.0) cf->SmoothObj(0.0);          // "Un-smooth" it
      cf->IntensityScaleObj(objmean/100.0);                                // "Un-scale" it
      cf->SaveRobj(clp->ObjOutFname());                                  
    }
  }
  catch(const std::exception& error) {
    cerr << "Error occured while writing output from fnirt" << endl;
    cerr << "Exception thrown with message: " << error.what() << endl; 
    return(EXIT_FAILURE);
  }

  freeparser(argc, argv);
  return(EXIT_SUCCESS);  // R.I.P.
}
    
