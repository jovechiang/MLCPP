#ifndef KMEANS_H
#define KMEANS_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <vector>
#include <set>

namespace KMeans{

	using namespace boost::numeric::ublas;

	typedef matrix<double> DataSet;
	typedef vector<double> Point;
	typedef unsigned int PointId;
	typedef unsigned int ClusterId;

	typedef std::set<PointId> Cluster;
	typedef std::vector<Cluster> Clusters;
	typedef std::vector<ClusterId> PointAssignment;
	typedef std::vector<Point> Centers;

	class kmeans{

		public: 
			void Initialization(int k);

			void KMeansClustering();
			
		private:
			void compute_centers();
			DataSet _dataset;
			Centers _centers;
			Clusters _clusters;
			PointAssignment _assign;

			
	};
};







#endif
