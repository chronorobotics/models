#ifndef MOMENT_ESTIMATOR_H
#define MOMENT_ESTIMATOR_H

#include <vector>

class CMoments;

class MomentEstimator {
	public:
		MomentEstimator(CMoments* parent_);
		~MomentEstimator();
		std::vector<double> sum_re;
		std::vector<double> sum_im;
		int count;

		void add_point(double phase);

		CMoments* get_parent() const;

	private:
		CMoments* parent;
};

#endif // MOMENT_ESTIMATOR_H
