/*
 *  first_newmat_vec.h
 *  
 *
 *  Created by Brian Patenaude on 12/08/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef NEWMAT_VEC_H
#define NEWMAT_VEC_H

#include <vector>
#include <newmat.h>

namespace FIRST_LIB{
	
	class first_newmat_vector{
		public:
			
	first_newmat_vector();
	~first_newmat_vector();
	
	
template<class T>
static std::vector<T> vectorToVector( const NEWMAT::Matrix & sm, const int & MaxModes);

template<class T>
static std::vector<T> vectorToVector( const NEWMAT::Matrix & sm);

template<class T>
static NEWMAT::ReturnMatrix vectorOfVectorsToMatrix( const std::vector< std::vector<T> > & vec);

template<class T>
static NEWMAT::ReturnMatrix vectorToDiagonalMatrix( const std::vector<T> & vec);

template<class T>
static std::vector< std::vector<T> > matrixToVector( const NEWMAT::Matrix & sm, const int & MaxModes);

template<class T>
static std::vector< std::vector<T> > matrixToVector( const NEWMAT::Matrix & sm);

template<class T>
static std::vector< std::vector<T> > matrixToVectorOfVectors(const NEWMAT::Matrix & m);

inline
static NEWMAT::ReturnMatrix unwrapMatrix(const NEWMAT::Matrix & m) 
{
	NEWMAT::ColumnVector munwrap(m.Nrows()*m.Ncols());
	unsigned int count=0;
	for (int i =0; i<m.Nrows() ; i++)
		for (int j =0; j<m.Ncols() ; j++,count++)
			munwrap.element(count)=m.element(i,j);
	return munwrap;
}

inline
static NEWMAT::ReturnMatrix wrapMatrix(const NEWMAT::ColumnVector & m) 
{
	
	NEWMAT::Matrix wrapped(m.Nrows()/3,3);
	unsigned int count=0;
	for (int i =0; i<m.Nrows() ; i+=3,count++)
	{
		wrapped.element(count,0)=m.element(i);
		wrapped.element(count,1)=m.element(i+1);
		wrapped.element(count,2)=m.element(i+2);
	}
	return wrapped;
}


	
	};
	
}

#endif
