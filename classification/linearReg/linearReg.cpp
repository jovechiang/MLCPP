#include "linearReg.h"

//featureSpace   m*n   matrix  ,  one row is one data entry , n dimension data
//labelVector    m*1   vector  

WeightVector linearRegDescGrad ( double alpha, unsigned int numIterations,
		const FeatureSpace & featureSpace, const LabelVector & labelVector)
{

	WeightVector weights(featureSpace.size2()) ;    //  n*1  vector	
	WeightVector h(featureSpace.size1());           // m*1 vector
	WeightVector err(featureSpace.size1());         //  m*1 vector

	//Initialize weights
	for(unsigned int i=0;i<weights.size(); ++i)
		weights(i) = 1.0;


	for(unsigned int i=0; i<numIterations; ++i){
		h = prod(featureSpace, weights);
		std::transform(h.begin, h.end(), h.begin(), sigmoid);

		err = h - labelVector;
		weights = weights + alpha * prod( trans(featureSpace) , err );
	}
	return weights;
}

WeightVector linearRegStocDescGrad(double alpha, unsigned int numIterations,
		const FeatureSpace & featureSpace, const LabelVector & labelVector)
{
	WeightVector weights(featureSpace.size2());      // n*1 vector
	double h;
	double e;


	//Initialize weights
	for(unsigned int i=0;i<weights.size(); ++i)
		weights(i) = 1.0;

	for(unsigned int k=0; k<numIterations; ++k){
		for(unsigned int i=0; i< FeatureSpace.size1(); ++i){
			matrix_row<matrix<double>> mi( featureSpace, i );
			h = inner_prod(mi, weights);
			h = sigmoid(h);
			e = labelVector(i) - h;
			weights = weights + alpha * e * mi;
		}
	}
	return weights;
}


