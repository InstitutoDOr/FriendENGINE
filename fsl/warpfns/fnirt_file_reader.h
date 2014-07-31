// Declarations of class used to decode and
// read files written by fnirt, and potentially
// by other pieces of software as long as they
// are valid displacement-field files.
//
// fnirt_file_reader.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2007 University of Oxford 
//

#ifndef fnirt_file_reader_h
#define fnirt_file_reader_h

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"

namespace NEWIMAGE {

class FnirtFileReaderException: public std::exception
{
private:
  std::string m_msg;
public:
  FnirtFileReaderException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("FnirtFileReader:: msg=" + m_msg).c_str();
  }

  ~FnirtFileReaderException() throw() {}
};

//
// Non-standard nifti intent-codes used for FSL-derived
// displacement/coefficient files.
//
#define FSL_FNIRT_DISPLACEMENT_FIELD       2006
#define FSL_CUBIC_SPLINE_COEFFICIENTS      2007
#define FSL_DCT_COEFFICIENTS               2008
#define FSL_QUADRATIC_SPLINE_COEFFICIENTS  2009

enum DisplacementFileType {FnirtSplineDispType, FnirtDCTDispType, FnirtFieldDispType, UnknownDispType, InvalidDispType};
enum AbsOrRelWarps {RelativeWarps, AbsoluteWarps, UnknownWarps};

// These are global routines that I have chosen to leave outside FnirtFileReader to
// allow them to be called independently. Theu probably don't really belong here, but
// I can't really come up with a better place for them.

void deffield2jacobian(const BASISFIELD::basisfield&   dx,
                       const BASISFIELD::basisfield&   dy,
                       const BASISFIELD::basisfield&   dz,
                       NEWIMAGE::volume<float>&        jac);

void deffield2jacobian(const BASISFIELD::basisfield&   dx,
                       const BASISFIELD::basisfield&   dy,
                       const BASISFIELD::basisfield&   dz,
                       const NEWMAT::Matrix&           aff,
                       volume<float>&                  jac);

void add_affine_part(const NEWMAT::Matrix&     aff,
                     unsigned int              indx,
                     NEWIMAGE::volume<float>&  warps);

void remove_affine_part(const NEWMAT::Matrix&     aff,
                        unsigned int              indx,
                        NEWIMAGE::volume<float>&  warps);

///////////////////////////////////////////////////////////////////////////////////////////////
//
// The FnirtFileReader is a class that reads, checks and decodes a displacement-field file.
// This can be a file created by fnirt, but it could also have been created for example
// from a fieldmap. For files created by other applications than fnirt it only reads 
// "proper" displacement fields. For fnirt-files it also reads coefficient-files.
//
// The ultimate purpose of FnirtFileReader and FnirtFileWriter is to provide an 
// interface that makes it easy and transparent not only to read and write displacement
// files, buat also to convert between the different representations. For example we
// want it to be easy to take a file with a direct representation of the field 
// (no basis-set) and write that out as a spline-coefficient file. Therefore we
// have internal space for all representations. When a file is read in it is stored
// in its "native" representation. When/if asked for another representation that
// representation is calculated and returned. 
//
///////////////////////////////////////////////////////////////////////////////////////////////

class FnirtFileReader
{
public:
  // Default constructor
  FnirtFileReader() : _type(InvalidDispType), _aor(UnknownWarps), _aff(IdentityMatrix(4)), _coef_rep(3) {}
  // Copy constructor
  FnirtFileReader(const FnirtFileReader& src);
  // Construct from file
  FnirtFileReader(const std::string& fname, AbsOrRelWarps wt=UnknownWarps, bool verbose=false)
  : _fname(fname)
  {
    common_read(fname,wt,verbose);
  }
  // Destructor. No clean up thanks to boost.
  ~FnirtFileReader() {}

  void Read(const std::string& fname, AbsOrRelWarps wt=UnknownWarps, bool verbose=false) {common_read(fname,wt,verbose);}
  bool IsValid() const {if (_type==InvalidDispType) return(false); else return(true);}
  std::string FileName() const {return(_fname);}
  DisplacementFileType Type() const {return(_type);}
  AbsOrRelWarps AbsOrRel() const {return(_aor);}
  std::vector<unsigned int> FieldSize() const;
  std::vector<double> VoxelSize() const;
  std::vector<unsigned int> KnotSpacing() const;
  unsigned int SplineOrder() const;
  std::vector<unsigned int> DCTOrder() const;
  NEWMAT::Matrix AffineMat() const {return(_aff);}
  NEWMAT::ReturnMatrix FieldAsNewmatMatrix(int indx=-1, bool inc_aff=false) const;
  NEWIMAGE::volume4D<float> FieldAsNewimageVolume4D(bool inc_aff=false) const;
  NEWIMAGE::volume<float> FieldAsNewimageVolume(unsigned int indx, bool inc_aff=false) const;
  NEWIMAGE::volume<float> Jacobian(bool inc_aff=false) const;
  BASISFIELD::splinefield FieldAsSplinefield(unsigned int indx, std::vector<unsigned int> ksp=std::vector<unsigned int>(0), unsigned int order=0) const;
  BASISFIELD::dctfield FieldAsDctfield(unsigned int indx, std::vector<unsigned int> order=std::vector<unsigned int>(0)) const;

protected:
  void common_read(const std::string& fname, AbsOrRelWarps wt, bool verbose);
  std::vector<boost::shared_ptr<BASISFIELD::basisfield> > read_coef_file(const NEWIMAGE::volume4D<float>&   vcoef,
                                                                         bool                               verbose) const;
  // void add_affine_part(NEWMAT::Matrix aff, unsigned int indx, NEWIMAGE::volume<float>& warps) const;

private:
  std::string                                                         _fname;
  DisplacementFileType                                                _type;
  AbsOrRelWarps                                                       _aor;
  NEWMAT::Matrix                                                      _aff;
  mutable std::vector<boost::shared_ptr<BASISFIELD::basisfield> >     _coef_rep;
  mutable boost::shared_ptr<NEWIMAGE::volume4D<float> >               _vol_rep;  
};

} // End namespace NEWIMAGE

#endif // end #ifndef fnirt_file_reader_h
