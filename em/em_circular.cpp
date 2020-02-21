#include <math.h>
#include "em_circular.h"

EMCircular::EMCircular() :
	cluster_count(),
	timestamps()
{

}

EMCircular::EMCircular(int cluster_count_) :
	cluster_count(cluster_count_),
	timestamps()
{

}

double EMCircular::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 604800.0f) / 604800;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

EMCircular::Timestamp::Timestamp(uint32_t time, int cluster_count, double weight_) :
	phase(time_to_phase(time)),
	weight(weight_),
	alpha()
{
	double s = 0;
	double r;

	for (int i = cluster_count; i; --i) {
		r = float(rand()) / RAND_MAX;
		s += r;
		alpha.push_back(r);
	}

	for (int i = cluster_count - 1; i >= 0; --i) {
		alpha[i] /= s;
	}
}

double EMCircular::get_loglikelihood() const {
	double result = 0;
	for (unsigned int i = 0; i < timestamps.size(); ++i) {
		result += log(get_density_at_d(timestamps[i].phase));
	}
	return result;
}

double EMCircular::get_density_at(uint32_t time) const {
	return get_density_at_d(time_to_phase(time));
}
