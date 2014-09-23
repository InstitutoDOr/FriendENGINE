//  
//  point_list.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2008 University of Oxford 
//

#ifndef point_list_h
#define point_list_h

#include <string>
#include <newmat.h>

namespace NEWIMAGE {

class PointListException: public std::exception
{
private:
  std::string m_msg;
public:
  PointListException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("PointList:: msg=" + m_msg).c_str();
  }

  ~PointListException() throw() {}
};

///////////////////////////////////////////////////////////////////////////////////////////////
//
// The PointList is a class that contains a list of points in a 3D-space. A PointList
// is constructed by reading a file that may be
//
// 1. A first output file
// 2. A text file with a coordinate triplet on each row/column.
// 
// The PointList will also be responsible for translating from the native formats of
// these files into "mm-coordinates", which is what we will ask PointList to output
// for us. Hence the constructor for each of those cases will also take an image
// file that defines the space within which the points are defined.
//
///////////////////////////////////////////////////////////////////////////////////////////////

class PointList
{
public:
  PointList(const std::string& pfname, const std::string& ifname);
  unsigned int NPoints() const {return(static_cast<unsigned int>(_points.Ncols()));}
  NEWMAT::ReturnMatrix Point(unsigned int indx, bool one_ext=false) const;
  NEWMAT::ReturnMatrix RawPoint(unsigned int indx, bool one_ext=false) const;
  NEWMAT::ReturnMatrix Vox2mm() const;
  NEWMAT::ReturnMatrix mm2Vox() const;
  void SetAffine(const NEWMAT::Matrix& aff);

private:
  int read_ascii_file(const std::string& fname);
  int read_first_file(const std::string& fname);
  bool is_maybe_vtk(const std::string& fname) const;

  std::string       _pfname;   // Name of file containing points
  std::string       _ifname;   // Name of image file defining space
  NEWMAT::Matrix    _SQform;   // qform/sform from image file
  NEWMAT::Matrix    _vox2mm;   // Voxel->mm-coordinates transform of _ifname
  NEWMAT::Matrix    _affine;   // Optional affine transform
  NEWMAT::Matrix    _points;   // Points in mm-coordinates of _ifname
}; 

} // End namespace NEWIMAGE

#endif // End #ifndef point_list_h
