#include <math.h>

#include "me_bicircular.h"

MEBicircular::Phase::Phase() :
	f1(),
	f2()
{

}

MEBicircular::Phase::Phase(uint32_t time) :
	f1(),
	f2()
{
	f1 = fmodf(time, 86400.0f) / 86400;
	if (f1 > 0.5) {
		f1 -= 1;
	}
	f1 *= M_PI * 2;

	f2 = fmodf(time, 604800.0f) / 604800;
	if (f2 > 0.5) {
		f2 -= 1;
	}
	f2 *= M_PI * 2;
}

MEBicircular::Phase::Phase(double f1_, double f2_) :
	f1(f1_),
	f2(f2_)
{

}

MEBicircular::Index::Index() :
	i1(),
	i2()
{

}

MEBicircular::Index::Index(int i1_, int i2_) :
	i1(i1_),
	i2(i2_)
{

}

MEBicircular::MEBicircular(int moment_count_) :
	moment_count(moment_count_),
	moment_indices()
{
	data.resize(moment_count*2, 0);

	int x = 1;
	int y = 0;
	for (int i = 0; i < moment_count; ++i) {
		moment_indices.push_back(Index(x, y));
		x--;
		y++;
		if (x < 0) {
			x = y;
			y = 0;
		}
	}
}

void MEBicircular::add_point(uint32_t time, double weight) {
	Phase phase(time);

	for (int i = 0; i < moment_count; ++i) {
		double foo = moment_indices[i].i1*phase.f1 + moment_indices[i].i2*phase.f2;
		data[2*i    ] += cos(foo) * weight;
		data[2*i + 1] += sin(foo) * weight;
	}

	count += weight;
}

std::vector<double> MEBicircular::get_moments() const {
	std::vector<double> result;
	result.resize(data.size());

	for (int i = data.size() - 1; i >= 0; --i) {
		result[i] = data[i] / count;
	}

	return result;
}

const std::vector<MEBicircular::Index>& MEBicircular::get_moment_indices() const {
	return moment_indices;
}

int MEBicircular::get_moment_count() const {
	return moment_count;
}
