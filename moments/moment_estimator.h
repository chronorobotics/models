#ifndef MOMENT_ESTIMATOR_H
#define MOMENT_ESTIMATOR_H

#include <vector>

class MomentEstimator {
	public:
		MomentEstimator();
		~MomentEstimator();
		std::vector<double> sum_re;
		std::vector<double> sum_im;
		int count;

		void add_point(double phase);
};

#endif // MOMENT_ESTIMATOR_H
