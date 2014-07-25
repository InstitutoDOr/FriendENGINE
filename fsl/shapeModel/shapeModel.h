/*
 *  shapeModel.h
 *  
 *
 *  Created by Brian Patenaude on 23/06/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
 #ifndef SHAPEMODEL_H
 #define SHAPEMODEL_H
 #include <vector>

namespace SHAPE_MODEL_NAME{


class shapeModel
{

public:
	
	std::vector< std::vector<unsigned int> > cells;
	std::vector< std::vector<unsigned int> > localTri;
	int NumberOfSubjects;
	shapeModel();
	shapeModel( const std::vector<float> & mshape, const std::vector< std::vector<float> > & modesshape, const std::vector<float> & se, \
				const std::vector<float> & ishape, const std::vector< std::vector<float> > & modesint, const std::vector<float> & ie, \
				const int & M, const std::vector<float> & errs);
				
	shapeModel( const std::vector<float> & mshape, const std::vector< std::vector<float> > & modesshape, const std::vector<float> & se, \
				const std::vector<float> & ishape, const std::vector< std::vector<float> > & modesint, const std::vector< std::vector<float> > & Iprec, const std::vector<float> & ie, \
				const int & M, const std::vector<float> & errs, const std::vector< std::vector<unsigned int> > & cellsin, const std::vector<int> & vlabels);
	
	shapeModel( const std::vector<float> &  mshape, const std::vector< std::vector<float> > & modesshape, const std::vector<float> & se, \
				const std::vector<float> & ishape, const std::vector< std::vector<float> > & modesint, const std::vector<float> & ie, \
				const int & M, const std::vector<float> & errs, const std::vector<short> & vmaskin);
	
	std::vector<float> getDeformedGrid( const std::vector<float>  & vars ) const ;
	
	std::vector<float> getDeformedIGrid( const std::vector<float> & vars) const ;
	 void printLabel(const unsigned int & i) const ;
	int getLabel( const unsigned int & i) const { return labels.at(i); }
	
	void registerModel(const std::vector< std::vector<float> > & flirtmat);
	std::vector<float> getOrigSpaceBvars(const std::vector<float> & bvars ) const;
	
	std::vector< std::vector<float> > registerModeVectors( const std::vector< std::vector<float> >& vmodes, const std::vector< std::vector<float> >& flirtmat);

	
	void setFoundMode(bool found) const { MODE_FOUND=found; }
	bool getFoundMode() const { return MODE_FOUND; } 
	void setMode(float m) const { mode=m; }
	float getMode() const { return mode; }
	bool getCondSet() const { return USE_COND; }
	void setCondSet( const bool & b ) const { USE_COND=b; }
	void setCondMats( const std::vector< std::vector<float> > & mat1, const std::vector< std::vector<float> > & mat2 , const unsigned int & k ){ kpred = k; condMat1 =mat1 ; condMat2=mat2; }

 std::vector< std::vector<float> > getCondMat1() const { return condMat1; }	
		std::vector< std::vector<float> > getCondMat2() const { return condMat2; }	
unsigned int getKPred() const { return kpred; }

	std::vector<float> smean;
	std::vector< std::vector<float> > smodes;
	std::vector< std::vector<float> > u_xfm;
	std::vector<float> d_xfm;
	
	std::vector< float > seigs;
	std::vector< float > sqrtseigs;

	std::vector< float > sqrtseigsi;

	std::vector< int > labels;

	std::vector<float> imean;
	std::vector< std::vector<float> > imodes;
	std::vector< float > ieigs;
	std::vector< std::vector<float> > i_precision;
	std::vector< float > Errs;
	std::vector<short> stmask;
	
	private:
		
		mutable std::vector< std::vector<float> > condMat1;
		mutable std::vector< std::vector<float> >  condMat2;
		mutable unsigned int kpred;
		mutable bool USE_COND;
	mutable bool STORE_REG_XFM;

		mutable bool MODE_FOUND;
		mutable float mode;
	
};

}

#endif

