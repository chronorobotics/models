#include "../CMoments.h"

#include "density_params.h"

#include "dp_von_mises.h"
#include "dp_double_von_mises.h"
#include "dp_sqdist.h"

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
		case DOUBLE_VON_MISES:
			return std::unique_ptr<DensityParams>(new DPDoubleVonMises(parent, count));
			break;
		case SQDIST:
			return std::unique_ptr<DensityParams>(new DPSqdist(parent, count));
			break;
		default:
			std::cerr << "Invalid Distribution " << (int) dist << std::endl;
			throw "DensityParams: Invalid distribution!";
			break;
	}
}

std::unique_ptr<DensityParams> DensityParams::load(CMoments* parent, FILE *file) {
	Distribution type;
	fread(&type, sizeof(Distribution), 1, file);
	std::unique_ptr<DensityParams> result = create(parent, type, 0);
	result->load(file);
	return std::move(result);
}

std::unique_ptr<DensityParams> DensityParams::importFromArray(CMoments *parent, double *array, int len, int &pos) {
	Distribution type = (Distribution) array[pos++];
	std::unique_ptr<DensityParams> result = create(parent, type, 0);
	result->importFromArray(array, len, pos);
	return std::move(result);
}
