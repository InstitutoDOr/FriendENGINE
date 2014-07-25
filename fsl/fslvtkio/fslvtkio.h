#ifndef _fslvtkIO
#define _fslvtkIO

#include <iostream>

#include <string>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include "meshclass/meshclass.h"
#include "first_lib/first_newmat_vec.h"

//using namespace std;
//using namespace mesh;


namespace fslvtkio{
		//using namespace mesh;
	
	class fslvtkIOException : public std::exception{
		
public:
		const char* errmesg;
		fslvtkIOException(const char* msg)
		{
			errmesg=msg;
		}
		
private:
			virtual const char* what() const throw()
		{
				cout<<errmesg<<endl;
				return errmesg;
		}
	};
	
	class fslvtkIO {
public:
		enum DataType{ POLYDATA, UNSTRUCTURED_GRID };
		fslvtkIO();
		fslvtkIO(const string & filename,const DataType i);
		
		~fslvtkIO();

		
		//----------------------SET COMMANDS-------------------------//
		void setPoints(const Matrix & pts);
		void setPoints(const vector<float> & pts);
		void setCells(const vector< vector<unsigned int> > & c){ Cells=c; }
		void setCells(const vector< vector<unsigned int> > & c, const vector<short> & c_t ){ Cells=c; Cell_Types=c_t; }
		
		void appendPointsAndPolygons(const Matrix & pts, const Matrix & Polygons);


		void setPolygons(const Matrix& m){ Polygons=m; }
		void setPolygons(const vector< vector<unsigned int> >& vm);
		void setMesh(const mesh::Mesh & m);
		//void setScalars(const Matrix & sc,const string & name);
		void setScalars(const Matrix& m){ Scalars=m; }
		
		void setVectors(const Matrix & vecs,const string &name);
		void setVectors(const Matrix& m){ Vectors=m; }

		template<class T>
		void setScalars(const vector<T> & sc);
		
			
		void addPointFieldData(const Matrix& m,const string & name, const string & type, const string & vtkAttType);			
		void addCellFieldData(const Matrix& m,const string & name, const string & type, const string & vtkAttType);
		
		template< class T > 
		void addFieldData(const vector<T> & m,const string & name, const string & type);
		
		void replaceFieldData(const Matrix& m,const string & name);

		void addFieldData(const Matrix& m,const string & name, const string & type);
		void addFieldData(vector< string >,string name);
		
		
		void printFieldDataNames() { for (vector<string>::iterator i=fieldDataNumName.begin(); i!=fieldDataNumName.end();i++) cout<<"field "<<*i<<endl; }

		
		//----------------------GET COMMANDS----------------------------//
		Matrix getPointsAsMatrix() const { return Points; }
		
		template<class T>
		std::vector<T> getPointsAsVector();
		
		template<class T>
		std::vector< std::vector<T> > getPointsAsVectorOfVectors(){ return FIRST_LIB::first_newmat_vector::matrixToVectorOfVectors<T>(Points); }
			
		template<class T>
		std::vector< std::vector<T> > getPolygonsAsVectorOfVectors(){ return FIRST_LIB::first_newmat_vector::matrixToVectorOfVectors<T>(Polygons); }

		ReturnMatrix getPolygons() const { return Polygons; };
		vector< vector<unsigned int> > getCells() const { return Cells; }
		vector<short> getCellTypes() const { return Cell_Types; }
								
		string getFieldName(const int & ind) const { return fieldDataNumName.at(ind); };
		unsigned int getNumberOfFields() const { return fieldDataNumName.size(); };
		Matrix getField(const string & name); 
		Matrix getField(const string & name, unsigned int & ind); 

		ReturnMatrix getScalars(){ return Scalars; }
		ReturnMatrix getVectors(){ return Vectors; }
		
		
		//----------------------I/O COMMANDS----------------------------//
		template<class T>
			void writePoints(ofstream & fshape, const string & str_typename);
		void writePolygons(ofstream & fshape);
		void writeCells(ofstream & fshape);
		//void writeUnstructuredGridCells(ofstream & fshape);
		void writeUnstructuredGridCellTypes(ofstream & fshape);
		void save(string s);
		
		void readUnstructuredGrid(string fname);
		void readPolyData(string fname);
		template<class T>
			void writePointData(ofstream & fshape, const string & str_typename );
		
		template <class T>
			void writeNumericField(ofstream & fvtk, const string & name, const string & type, const Matrix & Data);
		void writeStringField(ofstream & fvtk, const string & name, const vector<string> & v_string);
		
		void readFieldData(ifstream & fvtk);
		void readPointData(ifstream & fvtk, string & nextData);
		bool  readPoints(ifstream & fvtk);
		bool  readPolygons(ifstream & fvtk);
		
		template <class T>
			ReturnMatrix readField(ifstream & fvtk, const int & nrows,const int & mcols);
		
		void displayNumericFieldDataNames(); 				
				void displayNumericField(const string & name);

		//----------------------SETTING/ACCESS OF STATE VARIABLES----------------------------//
		void setSwitchRowsCols(bool n) { SWITCH_ROWS_COLS=n; }

		void setMAX(bool b) { MAX_SET=b; }
		void setMAX_Val( unsigned int n ) { MAX=n; }
		void setDataType(const DataType & dtype ){ dt=dtype; }
		void setBinaryWrite(bool state) { BINARY=state; }
		bool getBinaryWrite() { return BINARY; }
		static bool myExceptions(int e);
		
protected:
		Matrix Scalars;
		Matrix Vectors;
		Matrix Points;
		//cell data
		Matrix Polygons;
		
		
		
private:
#define BINFLAG 42
			//----------------------STATE VARIABLES/CONSTANTS----------------------------//
				
#ifdef PPC64
	int m_n;
#endif
			bool BINARY;//state variable for read write
			bool SWAP_BYTES;//only used for binary read
				bool MAX_SET; //is a max number of columns imposed
				bool SWITCH_ROWS_COLS;//SWITCHES THE COLUMN/ROWS , USED FOR FIXING OLD FILES

				unsigned int ST_COUNT;
	
					unsigned int MAX;
					DataType dt;//i.e. POLYDATA
						
						
						
						
						//----------------------READ COMMANDS----------------------------//
							
						//----------------------DATA----------------------------//
						
						//point data
						string scalarsName, vectorsName;
						
						vector< vector<unsigned int> > Cells;
						vector<short> Cell_Types;
						
						string cd_scalarsName, cd_vectorsName;
						Matrix cd_Scalars;
						Matrix cd_Vectors;
						
						
						//field data
						vector< Matrix > fieldDataNum;
						vector< string > fieldDataNumName;
						vector< string > fieldDataNumType;
						
						vector< vector<string> > fieldDataStr;
						vector< string > fieldDataStrName;
						
						//defaults FIELDData
						vector< string > pd_list;
						vector< string > pd_type;
						
						vector< string > cd_list;
						vector< string > cd_type;
						
	};
}
#endif
