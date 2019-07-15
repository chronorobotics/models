#ifndef MOMENT_ESTIMATOR_H
#define MOMENT_ESTIMATOR_H


class MomentEstimator {
	public:
		MomentEstimator();
		~MomentEstimator();
		double* sum_re;
		double* sum_im;
		int count;

		void add_point(double phase);
};

#endif // MOMENT_ESTIMATOR_H
