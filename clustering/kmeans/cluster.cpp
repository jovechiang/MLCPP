#include "cluster.h"

namespace KMeans{

	void kmeans::Initialization(int k){
		_clusters.resize(k);
		_assign.resize(_dataset.size1());
		_centers.resize(k);
		ClusterId cid;

		for(PointId pid = 0; pid < _dataset.size1() ; ++pid){
			cid = pid % k;
			_clusters[cid].insert(pid);
			_assign[pid] = cid;
		}
	}

	void KMeansClustering(){
		bool move = true;
		while(move){
			move = false;
			compute_centers();

			for(PointId pid; pid < _dataset.size1(); ++pid){
				matrix_row<matrix<double> > mi(_dataset, _clusters[cid][pid]);
				Point diff = trans(mi) - _centers[_assign[pid]];
				double mindistance = norm_2(diff);

				//find the closest center
				ClusterId closest = _assign[pid];
				for(ClusterId cid; cid < _centers.size(); ++cid){
					if(cid == _assign[pid]){
						continue;
					}

					double thisDistance = norm_2(trans(mi) - _centers[cid]);
					if(mindistance > thisDistance ){
						mindistance = thisDistance;
						move = true;
						closest = cid;
					}

					if(move){
						_clusters[_assign[pid]].erase(pid);
						_assign[pid] = closest;
						_clusters[closest].insert(pid);
					}

				}
			}

		}
	}

	void compute_centers(){
		
		for(ClusterId cid = 0; cid < _clusters.size(); ++cid){
			Point p(_dataset.size2());
			for(unsigned int i=0; i< p.size(); ++i){
				p(i) = 0.0;
			}

			for(PointId pid = 0 ; pid < _clusters[cid].size(); ++pid){
				matrix_row<matrix<double> > mi(_dataset, _clusters[cid][pid]);
				p += trans(mi);
			}
			p/=_clusters[cid].size();
			_centers[cid] = p;
		}
	}
	
};
