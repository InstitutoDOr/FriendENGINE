/*
 *  first_newmat_vec.cpp
 *  
 *
 *  Created by Brian Patenaude on 12/08/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "first_newmat_vec.h"

using namespace std;
using namespace NEWMAT;

namespace FIRST_LIB{

first_newmat_vector::first_newmat_vector(){

}

first_newmat_vector::~first_newmat_vector(){

}

template<class T>
std::vector<T> first_newmat_vector::vectorToVector( const Matrix & sm, const int & MaxModes)
{
	std::vector<T> vecM;	
	if (sm.Nrows()==1){
		for (int i=0;i<  MaxModes ; i++){
			vecM.push_back(static_cast<T>(sm.element(0,i)));
		}
	}else{
		for (int i=0;i<  MaxModes ; i++){
			vecM.push_back(static_cast<T>(sm.element(i,0)));
		}
		
		
	}
	
	return vecM;
}

template std::vector<unsigned int> first_newmat_vector::vectorToVector<unsigned int>( const Matrix & sm, const int & MaxModes);
template std::vector<int> first_newmat_vector::vectorToVector<int>( const Matrix & sm, const int & MaxModes);
template std::vector<float> first_newmat_vector::vectorToVector<float>( const Matrix & sm, const int & MaxModes);
template std::vector<double> first_newmat_vector::vectorToVector<double>( const Matrix & sm, const int & MaxModes);

template<class T>
std::vector<T> first_newmat_vector::vectorToVector( const Matrix & sm)
{
	vector<T> vecM;	
	if (sm.Nrows()==1){
		for (int i=0;i< sm.Ncols() ; i++){
			vecM.push_back(static_cast<T>(sm.element(0,i)));
		}
	}else{
		for (int i=0;i< sm.Nrows() ; i++){
			vecM.push_back(static_cast<T>(sm.element(i,0)));
		}
		
		
	}
	
	return vecM;
}

template std::vector<unsigned int> first_newmat_vector::vectorToVector<unsigned int>( const Matrix & sm);
template std::vector<int> first_newmat_vector::vectorToVector<int>( const Matrix & sm);
template std::vector<float> first_newmat_vector::vectorToVector<float>( const Matrix & sm);
template std::vector<double> first_newmat_vector::vectorToVector<double>( const Matrix & sm);


template<class T>
ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix( const vector< vector<T> > & vec)
{
	Matrix out(vec.size(), vec.at(0).size());
	unsigned int row=0;

	for (typename vector< vector<T> >::const_iterator i=vec.begin() ; i!=vec.end();i++, row++)
	{
		unsigned int col=0;
		for (typename vector<T>::const_iterator j=i->begin() ; j!=i->end();j++,col++)
			out.element(row,col)=*j;
	}		
	out.Release();
	return out;
}

template ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix<unsigned int>( const vector< vector<unsigned int> > & vec);
template ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix<int>( const vector< vector<int> > & vec);
template ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix<float>( const vector< vector<float> > & vec);
template ReturnMatrix first_newmat_vector::vectorOfVectorsToMatrix<double>( const vector< vector<double> > & vec);

template<class T>
ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix( const vector<T> & vec)
{
	DiagonalMatrix out(vec.size());
	unsigned int row=0;

	for (typename vector<T> ::const_iterator i=vec.begin() ; i!=vec.end();i++, row++)
		out.element(row)=*i;
		
	out.Release();
	return out;
}

template ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix<unsigned int>( const vector<unsigned int> & vec);
template ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix<int>( const vector<int> & vec);
template ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix<float>( const vector<float> & vec);
template ReturnMatrix first_newmat_vector::vectorToDiagonalMatrix<double>( const vector<double> & vec);


template<class T>
vector< vector<T> > first_newmat_vector::matrixToVector( const Matrix & sm, const int & MaxModes)
{
	vector< vector<T> > vecM;	
    for (int j=0;j< MaxModes ; j++){
		vector<T> mode;
	    for (int i=0;i< sm.Nrows() ; i++){
			mode.push_back(sm.element(i,j));
		}
		vecM.push_back(mode);
	}   
	
	return vecM;
}

template vector< vector<unsigned int> > first_newmat_vector::matrixToVector<unsigned int>( const Matrix & sm, const int & MaxModes);
template vector< vector<int> > first_newmat_vector::matrixToVector<int>( const Matrix & sm, const int & MaxModes);
template vector< vector<float> > first_newmat_vector::matrixToVector<float>( const Matrix & sm, const int & MaxModes);
template vector< vector<double> > first_newmat_vector::matrixToVector<double>( const Matrix & sm, const int & MaxModes);


template<class T>
vector< vector<T> > first_newmat_vector::matrixToVector( const Matrix & sm)
{
	vector< vector<T> > vecM;	
    for (int j=0;j< sm.Ncols() ; j++)
	{
		vector<T> mode;
	    for (int i=0;i< sm.Nrows() ; i++)
			mode.push_back(static_cast<T>(sm.element(i,j)));
		vecM.push_back(mode);
	}   

	return vecM;
}

template vector< vector<unsigned int> > first_newmat_vector::matrixToVector<unsigned int>( const Matrix & sm);
template vector< vector<int> > first_newmat_vector::matrixToVector<int>( const Matrix & sm);
template vector< vector<float> > first_newmat_vector::matrixToVector<float>( const Matrix & sm);
template vector< vector<double> > first_newmat_vector::matrixToVector<double>( const Matrix & sm);

template<class T>
std::vector< std::vector<T> > first_newmat_vector::matrixToVectorOfVectors(const NEWMAT::Matrix & m){
	vector< vector<T> > all;
	for (int i=0;i<m.Nrows();i++)
	{
		vector<T> row;
		for (int j=0; j<m.Ncols();j++)
			row.push_back(static_cast<T>(m.element(i,j)));
		all.push_back(row);
	}
	return all;
}
template std::vector< std::vector<short> > first_newmat_vector::matrixToVectorOfVectors<short>(const NEWMAT::Matrix & m);
template std::vector< std::vector<unsigned int> > first_newmat_vector::matrixToVectorOfVectors<unsigned int>(const NEWMAT::Matrix & m);
template std::vector< std::vector<int> > first_newmat_vector::matrixToVectorOfVectors<int>(const NEWMAT::Matrix & m);
template std::vector< std::vector<float> > first_newmat_vector::matrixToVectorOfVectors<float>(const NEWMAT::Matrix & m);
template std::vector< std::vector<double> > first_newmat_vector::matrixToVectorOfVectors<double>(const NEWMAT::Matrix & m);




}
