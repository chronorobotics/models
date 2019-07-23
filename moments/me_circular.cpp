#include <math.h>

#include "me_circular.h"

MECircular::MECircular(int moment_count_) :
	MomentEstimator(),
	moment_count(moment_count_)
{
	data.resize(moment_count*2);
}

MECircular::~MECircular()
{

}

void MECircular::add_point(uint32_t time) {
	double phase = time_to_phase(time);

	for (int i = moment_count - 1; i >= 0; --i) {
		data[2*i    ] += cos((i+1)*phase);
		data[2*i + 1] += sin((i+1)*phase);
	}
	count++;
}

double MECircular::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 86400.0f) / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

std::vector<double> MECircular::get_moments() const {
	std::vector<double> result;
	result.resize(data.size());

	for (int i = data.size() - 1; i >= 0; --i) {
		result[i] = data[i] / count;
	}

	return result;
}

int MECircular::get_moment_count() const {
	return data.size()/2;
}
