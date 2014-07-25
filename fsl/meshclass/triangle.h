/*  Copyright (C) 1999-2004 University of Oxford  */

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

#ifndef _triangle
#define _triangle
#include <iostream>
#include <vector>

#include"point.h"

using namespace std;

namespace mesh{

class Mpoint;

class Triangle
{
 private:
  Mpoint * _vertice[3];
  float _value;
  
 public:
  Triangle(Mpoint* const p1,Mpoint* const p2,Mpoint* const p3, float val=0);
  ~Triangle(void);
  Triangle(Triangle & t);          //prevents from using copy constructor
  Triangle operator=(Triangle &t); //prevents from using affectation operator

  vector<double> data; //can be used to store extra-data attached to the triangle

  const Pt centroid() const;
  const Vec normal() const;
  const Vec area(const Mpoint* p) const; 
  Mpoint * const get_vertice(const int i) const;
  const float get_value() const {return _value;};
  void set_value(const float val){_value=val;};
  void swap();             //changes triangle orientation
  const int operator <(const Triangle *const t) const;   //checks if two triangles are adjacents and well-oriented
  const bool operator ==(const Triangle & t) const;
  const bool intersect(const Triangle & t) const; //checks if two triangles intersect

  bool oriented;           //has the triangle been well oriented ?
  

};

}
#endif 
