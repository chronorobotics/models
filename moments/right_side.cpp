#include "../CMoments.h"
#include "moment_estimator.h"

#include "right_side.h"

RightSide::RightSide() :
	moment_re(),
	moment_im(),
	count(CMoments::moment_count)
{
	moment_re = new double[CMoments::moment_count];
	moment_im = new double[CMoments::moment_count];
}

RightSide::RightSide(const MomentEstimator& me) :
	moment_re(),
	moment_im(),
	count(CMoments::moment_count)
{
	moment_re = new double[CMoments::moment_count];
	moment_im = new double[CMoments::moment_count];

	for (int i = 0; i < CMoments::moment_count; ++i) {
		moment_re[i] = me.sum_re[i] / me.count;
		moment_im[i] = me.sum_im[i] / me.count;
	}
}

RightSide::~RightSide() {
	delete[] moment_re;
	delete[] moment_im;
}
