// Declarations of class used for point-matching in FNIRT
//
// matching_points.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2008 University of Oxford 
//


#ifndef matching_points_h
#define matching_points_h

#include <string>
#include "newmat.h"
#include "miscmaths/bfmatrix.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "warpfns/point_list.h"

namespace FNIRT {

class MatchingPointsException: public std::exception
{
private:
  std::string m_msg;
public:
  MatchingPointsException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return(string("MatchingPoints:: " + m_msg).c_str());
  }

  ~MatchingPointsException() throw() {}
};

class MatchingPoints
{
public:
  MatchingPoints(const NEWIMAGE::PointList&  points_in_ref, const NEWIMAGE::PointList& points_in_in) 
  {
    if (points_in_ref.NPoints() != points_in_in.NPoints()) throw MatchingPointsException("MatchingPoints: Size mismatch between point-lists");
    _in_ref = boost::shared_ptr<NEWIMAGE::PointList>(new NEWIMAGE::PointList(points_in_ref));
    _in_in = boost::shared_ptr<NEWIMAGE::PointList>(new NEWIMAGE::PointList(points_in_in));
  }
  virtual ~MatchingPoints() {}

  // Returns the difference in position in dimension dim for point-pair given by index
  virtual double Diff(unsigned int indx, unsigned int dim, const BASISFIELD::basisfield& bf) const;
  // Returns the sum of squared differences in position along dimension dim
  virtual double SSD(unsigned int dim, const BASISFIELD::basisfield& bf) const;
  // Returns the gradient of the sum-of-squared differences for the dimension given by dim
  virtual NEWMAT::ReturnMatrix SSD_Gradient(unsigned int dim, const BASISFIELD::basisfield& bf) const;
  // Returns the Hessian of the sum-of-squared differences for the dimension given by dim
  virtual boost::shared_ptr<MISCMATHS::BFMatrix> SSD_Hessian(unsigned int                     dim, 
                                                             const BASISFIELD::basisfield&    bf,
                                                             MISCMATHS::BFMatrixPrecisionType prec=BFMatrixDoublePrecision) const;
 
private:
  boost::shared_ptr<NEWIMAGE::PointList>  _in_ref;
  boost::shared_ptr<NEWIMAGE::PointList>  _in_in;
};

}

#endif // end #ifndef matching_points_h
