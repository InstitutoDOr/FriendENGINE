/* first_mesh.cpp
Created by Brian Patenaude on 12/08/2008.
*/
#ifndef FIRST_MESH_H
#define FIRST_MESH_H

#include <vector>
#include "newimage/newimageall.h"

namespace FIRST_LIB{

class first_mesh{

	public:
		first_mesh();
		~first_mesh();
static void getBounds(const std::vector<float> & m, int *bounds, const float & xdim, const float & ydim, const float & zdim);
static void draw_segment(NEWIMAGE::volume<short>& image, const float & p1x,const float & p1y,const float & p1z,  const float & p2x,const float & p2y,const float & p2z, int label);
static NEWIMAGE::volume<short> draw_mesh(const NEWIMAGE::volume<short> & image, const std::vector<float> & pts, const std::vector< std::vector<unsigned int> > & triangles, int label);
static NEWIMAGE::volume<short> make_mask_from_mesh(const NEWIMAGE::volume<float> & image, const std::vector<float> & m, const std::vector< std::vector<unsigned int> > & triangles, const int & label, const int * bounds);


template<class T>
static vector< vector<T> >  findNeighbours(const vector< vector<T> > & cells, const unsigned int & N);

template<class T>
static vector< vector<T> > findNeighbourTriangles(const vector< vector<T> > & cells, const unsigned int & N);



template<class T, class T2>
static void maxTriangle(const vector<T> & pts, const vector< vector<T2> > & localTri, const vector< vector<T2> > & cells,vector<T> & vx, vector<T> & vy, vector<T> & vz );

template<class T,class T2>
static bool self_intersection_test(const vector< vector<T2> > & tri_cells, const vector<T> & pts);

/*template<class T>
short pointTest(const T & px_0, const T & py_0, const T & pz_0, const T & px_1, const T & py_1, const T & pz_1, \
				const T & px_2, const T & py_2, const T & pz_2, const T & nx, const T & ny, const T & nz. \
				T & p0_x, T & p0_y, T & p0_z, T & dx, T & dy, T & dz);*/
template<class T> 
static bool pt_equal(const T* pt1, const T* pt2)
{
	return ((pt1[0]==pt2[0]) && (pt1[1]==pt2[1]) && (pt1[2]==pt2[2]) ) ;
}

template<class T> 
static T* vec_sub(const T* vec2, const T* vec1)
{
	T* vsub = new T[3];
	vsub[0]=vec2[0]-vec1[0];
	vsub[1]=vec2[1]-vec1[1];
	vsub[2]=vec2[2]-vec1[2];
	
	return vsub;
}


template<class T> 
static T dot_prod(const T* vec1, const T* vec2 )
{	 
	return ( vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2] );
}

template<class T> 
static T* cross_prod(const T* v1, const T* v2)
{
//	cout<<"cross "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<" "<<v2[0]<<" "<<v2[1]<<" "<<v2[2]<<endl;
	T* norm = new T[3];
	norm[0]=(v1[1]*v2[2]-v1[2]*v2[1]);
	norm[1]=(v1[2]*v2[0]-v1[0]*v2[2]);
	norm[2]=(v1[0]*v2[1]-v1[1]*v2[0]);

	return norm;
}

template<class T> 
static bool triangle_intersection(const T* V10, const T* V11, const T* V12, const T* V20, const T* V21, const T* V22 );

template<class T,class T2>
static void normal(const vector<T> & pts, const vector< vector<T2> > & localTri, const vector< vector<T2> > & cells,vector<T> & vx, vector<T> & vy, vector<T> & vz );


template<class T,class T2>
static void  medium_neighbours(const vector<T> & pts, const vector< vector<T2> > & neighbours, const vector< vector<T2> > & cells,\
									vector<T> & vx, vector<T> & vy, vector<T> & vz );						  
};

}
#endif
