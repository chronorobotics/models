#include "../CMoments.h"
#include "moment_estimator.h"

#include "right_side.h"

RightSide::RightSide(CMoments* parent_) :
	moment_re(),
	moment_im(),
	count(parent_->get_moment_count()),
	parent(parent_)
{
	moment_re.resize(count);
	moment_im.resize(count);
}

RightSide::RightSide(const MomentEstimator& me) :
	moment_re(),
	moment_im(),
	count(me.get_parent()->get_moment_count()),
	parent(me.get_parent())
{
	moment_re.resize(count);
	moment_im.resize(count);

	for (int i = 0; i < count; ++i) {
		moment_re[i] = me.sum_re[i] / me.count;
		moment_im[i] = me.sum_im[i] / me.count;
	}
}

RightSide::~RightSide() {

}
