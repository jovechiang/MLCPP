#ifndef LINEARREG_H
#define LINEARREG_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

typedef matrix<double> FeatureSpace;
typedef vector<unsigned short> LabelVector;
typedef vector<double> WeightVector;

inline double sigmoid(double x) {return 1.0/ (1 + exp( -x )); }

WeightVector linearRegDescGrad(double alpha, unsigned int numIterations,
	   const FeatureSpace & featureSpace, const LabelVector & labelVector);

WeightVector linearRegStocDescGrad(double alpha, unsigned int numIterations,
		const FeatureSpace & featureSpace, const LabelVector & labelVector);

#endif 
