#include <math.h>
#include "em_circular.h"

EMCircular::EMCircular(int cluster_count_, double period_) :
	cluster_count(cluster_count_),
	period(period_),
	timestamps()
{

}

double EMCircular::time_to_phase(uint32_t time, double period_) {
	float phase = fmodf(time, period_) / period_;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

EMCircular::Timestamp::Timestamp(uint32_t time, int cluster_count, double period_) :
	phase(time_to_phase(time, period_)),
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
