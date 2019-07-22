#include "../CMoments.h"

#include "density_params.h"

DensityParams::DensityParams(CMoments* parent_, int count_) :
	count(),
	parent(parent_)
{
	if (count_ < 0) {
		count = parent->get_cluster_count();
	} else {
		count = count_;
	}
}

DensityParams::~DensityParams() {

}
