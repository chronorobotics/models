#include "math.h"
#include "../CMoments.h"

#include "moment_estimator.h"

MomentEstimator::MomentEstimator(CMoments* parent_) :
	sum_re(),
	sum_im(),
	count(0),
	parent(parent_)
{
	sum_re.resize(parent->get_moment_count(), 0);
	sum_im.resize(parent->get_moment_count(), 0);
}

MomentEstimator::~MomentEstimator() {

}

void MomentEstimator::add_point(double phase) {
	for (int i = parent->get_moment_count()-1; i >= 0; --i) {
		sum_re[i] += cos((i+1)*phase);
		sum_im[i] += sin((i+1)*phase);
	}
	count++;
}

CMoments* MomentEstimator::get_parent() const {
	return parent;
}
