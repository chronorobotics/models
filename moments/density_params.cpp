#include "../CMoments.h"

#include "density_params.h"

DensityParams::DensityParams(int count_) :
	count()
{
	if (count_ < 0) {
		count = CMoments::moment_count*2/3;
	} else {
		count = count_;
	}
}

DensityParams::~DensityParams() {

}
