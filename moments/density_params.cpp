#include "../CMoments.h"

#include "density_params.h"

#include "dp_von_mises.h"

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

std::unique_ptr<DensityParams> DensityParams::create(CMoments* parent, Distribution dist, int count) {
	switch (dist) {
		case VON_MISES:
			return std::unique_ptr<DensityParams>(new DPVonMises(parent, count));
			break;
		default:
			throw "DensityParams: Invalid distribution!";
			break;
	}
}
