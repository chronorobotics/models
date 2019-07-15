#include "math.h"
#include "../CMoments.h"

#include "moment_estimator.h"

MomentEstimator::MomentEstimator() :
	sum_re(),
	sum_im(),
	count(0)
{
	sum_re = new double[CMoments::moment_count]();
	sum_im = new double[CMoments::moment_count]();
}

MomentEstimator::~MomentEstimator() {
	delete[] sum_re;
	delete[] sum_im;
}

void MomentEstimator::add_point(double phase) {
	for (int i = 1; i <= CMoments::moment_count; ++i) {
		sum_re[i] += cos(i*phase);
		sum_im[i] += sin(i*phase);
	}
	count++;
}
