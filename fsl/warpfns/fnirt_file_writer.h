// Declarations of class used to write fnirt
// displacement files.
//
// fnirt_file_writer.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2007 University of Oxford 
//

#ifndef fnirt_file_writer_h
#define fnirt_file_writer_h

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "basisfield/dctfield.h"
#include "fnirt_file_reader.h"

namespace NEWIMAGE {

class FnirtFileWriterException: public std::exception
{
private:
  std::string m_msg;
public:
  FnirtFileWriterException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("FnirtFileReader:: msg=" + m_msg).c_str();
  }

  ~FnirtFileWriterException() throw() {}
};

class FnirtFileWriter
{
public:
  // Constructor for coefficient file
  FnirtFileWriter(const std::string&                                               fname, 
                  const std::vector<boost::shared_ptr<BASISFIELD::basisfield> >&   fields,
                  NEWMAT::Matrix                                                   aff=NEWMAT::IdentityMatrix(4));
  // "Simpler" constructor for coefficient file
  FnirtFileWriter(const std::string&                   fname, 
                  const BASISFIELD::basisfield&        fieldx,
                  const BASISFIELD::basisfield&        fieldy,
                  const BASISFIELD::basisfield&        fieldz,
                  NEWMAT::Matrix                       aff=NEWMAT::IdentityMatrix(4))
  {
    common_coef_construction(fname,fieldx,fieldy,fieldz,aff);
  }
  // Constructor for field file
  FnirtFileWriter(const std::string&                   fname,
                  const NEWIMAGE::volume<float>&       ref,
                  const NEWIMAGE::volume4D<float>&     vfields,                         
                  NEWMAT::Matrix                       aff=NEWMAT::IdentityMatrix(4));
  // Another constructor for field file
  FnirtFileWriter(const std::string&                     fname,
                  const NEWIMAGE::volume4D<float>&       vfields);
  // Another constructor for field file
  FnirtFileWriter(const std::string&                   fname,
		  const NEWIMAGE::volume<float>&       fieldx,                         
		  const NEWIMAGE::volume<float>&       fieldy,                         
		  const NEWIMAGE::volume<float>&       fieldz);
  // Another constructor for field file
  FnirtFileWriter(const std::string&                   fname,
                  const NEWIMAGE::volume<float>&       ref,
                  const NEWIMAGE::volume<float>&       fieldx,                         
                  const NEWIMAGE::volume<float>&       fieldy,                         
                  const NEWIMAGE::volume<float>&       fieldz,                         
                  NEWMAT::Matrix                       aff=NEWMAT::IdentityMatrix(4))
  {
    common_field_construction(fname,ref,fieldx,fieldy,fieldz,aff);
  }

protected:
  void common_coef_construction(const std::string&                   fname, 
                                const BASISFIELD::basisfield&        fieldx,
                                const BASISFIELD::basisfield&        fieldy,
                                const BASISFIELD::basisfield&        fieldz,
                                const NEWMAT::Matrix&                aff);

  void common_field_construction(const std::string&                 fname,
                                 const NEWIMAGE::volume<float>&     ref,
                                 const NEWIMAGE::volume<float>&     fieldx,
                                 const NEWIMAGE::volume<float>&     fieldy,
                                 const NEWIMAGE::volume<float>&     fieldz, 
				 const NEWMAT::Matrix&              aff);
private:
};

} // End namespace NEWIMAGE

#endif // End #ifndef fnirt_file_writer_h
