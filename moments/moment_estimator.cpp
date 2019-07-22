#include "math.h"
#include "../CMoments.h"

#include "moment_estimator.h"

MomentEstimator::MomentEstimator() :
	sum_re(),
	sum_im(),
	count(0)
{
	sum_re.resize(CMoments::moment_count, 0);
	sum_im.resize(CMoments::moment_count, 0);
}

MomentEstimator::~MomentEstimator() {

}

void MomentEstimator::add_point(double phase) {
	for (int i = 0; i < CMoments::moment_count; ++i) {
		sum_re[i] += cos((i+1)*phase);
		sum_im[i] += sin((i+1)*phase);
	}
	count++;
}
