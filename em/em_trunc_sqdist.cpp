#include "em_trunc_sqdist.h"

#include <math.h>
#include <iostream>
#include <sstream>

EMTruncSqdist::EMTruncSqdist(int cluster_count_) :
	cluster_count(cluster_count_),
	clusters(),
	timestamps()
{
	double s = 0;
	for (int i = 0; i < cluster_count; ++i) {
		clusters.push_back(Cluster());
		s += clusters[i].weight;
		double r = double(rand())/RAND_MAX;
		double f = 2*M_PI*double(rand())/RAND_MAX;
		clusters[i].xx = r*cos(f);
		clusters[i].yy = r*sin(f);
		clusters[i].theta = M_PI*double(rand())/RAND_MAX;
		clusters[i].precalc();
	}

	for (int i = 0; i < cluster_count; ++i) {
		clusters[i].weight /= s;
	}
}

void EMTruncSqdist::expectation() {
	std::cerr << "Performing expectation ..." << std::endl;

	for (int i = 0; i < timestamps.size(); ++i) {
		double s = 0;
		for (int j = 0; j < clusters.size(); ++j) {
			double a = clusters[j].density_at(timestamps[i].phase);
			timestamps[i].alpha[j] = a;
			s += a;
		}

		for (int j = 0; j < clusters.size(); ++j) {
			timestamps[i].alpha[j] /= s;
		}
	}
}

double EMTruncSqdist::maximisation() {
	double shift = 0;
	std::cerr << "Performing maximisation ..." << std::endl;

	for (int i = 0; i < clusters.size(); ++i) {
		double last_xx = clusters[i].xx;
		double last_yy = clusters[i].yy;
		double last_theta = clusters[i].theta;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			s += timestamps[j].alpha[i];
		}
		clusters[i].weight = s / timestamps.size();

		double m1r = 0;
		double m1i = 0;
		double m2r = 0;
		double m2i = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			m1r += cos(timestamps[j].phase) * timestamps[j].alpha[i];
			m1i += sin(timestamps[j].phase) * timestamps[j].alpha[i];
			m2r += cos(timestamps[j].phase*2) * timestamps[j].alpha[i];
			m2i += sin(timestamps[j].phase*2) * timestamps[j].alpha[i];
		}
		m1r /= timestamps.size();
		m1i /= timestamps.size();
		m2r /= timestamps.size();
		m2i /= timestamps.size();
		clusters[i].estimate(m1r, m1i, m2r, m2i);


		double delta_xx = last_xx - clusters[i].xx;
		double delta_yy = last_yy - clusters[i].yy;
		double delta_theta = last_theta - clusters[i].theta;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_xx*delta_xx + delta_yy*delta_yy + delta_theta*delta_theta + delta_weight*delta_weight;
	}

	return sqrt(shift);
}

double EMTruncSqdist::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 604800.0f) / 604800;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

void EMTruncSqdist::train() {
	double shift = 555;
	do {
		expectation();
		shift = maximisation();
		std::cout << "shift = " << shift << std::endl;
	} while (isnan(shift) || shift > 1E-2);
}

void EMTruncSqdist::add_time(uint32_t time) {
	timestamps.push_back(Timestamp(time, cluster_count));
}

EMTruncSqdist::Timestamp::Timestamp(uint32_t time, int cluster_count) :
	phase(time_to_phase(time)),
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

EMTruncSqdist::Cluster::Cluster() :
	xx(),
	yy(),
	theta(),
	weight(float(rand()) / RAND_MAX),
	corrective(),
	trunc()
{
	double r = double(rand())/RAND_MAX;
	double f = 2*M_PI*double(rand())/RAND_MAX;
	xx = r*cos(f);
	yy = r*sin(f);
	theta = M_PI*double(rand())/RAND_MAX;
}

EMTruncSqdist::Cluster::Cluster(double xx_, double yy_, double theta_, double weight_) :
	xx(xx_),
	yy(yy_),
	theta(theta_),
	weight(),
	corrective(),
	trunc()
{
	precalc();
}

void EMTruncSqdist::Cluster::print() {
	std::cout << "(" << xx << ", " << yy << ", " << theta << ", " << weight << ")";
}

void EMTruncSqdist::Cluster::save(FILE* file, bool lossy) {
	fwrite(&xx, sizeof(double), 1, file);
	fwrite(&yy, sizeof(double), 1, file);
	fwrite(&theta, sizeof(double), 1, file);
	fwrite(&weight, sizeof(double), 1, file);
	fwrite(&corrective, sizeof(double), 1, file);
	fwrite(&trunc, sizeof(double), 1, file);
}

void EMTruncSqdist::Cluster::load(FILE* file) {
	fread(&xx, sizeof(double), 1, file);
	fread(&yy, sizeof(double), 1, file);
	fread(&theta, sizeof(double), 1, file);
	fread(&weight, sizeof(double), 1, file);
	fread(&corrective, sizeof(double), 1, file);
	fread(&trunc, sizeof(double), 1, file);
}

void EMTruncSqdist::Cluster::importFromArray(double* array, int len, int& pos) {
	xx = array[pos++];
	yy = array[pos++];
	theta = array[pos++];
	weight = array[pos++];
	corrective = array[pos++];
	trunc = array[pos++];
}

void EMTruncSqdist::Cluster::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = xx;
	array[pos++] = yy;
	array[pos++] = theta;
	array[pos++] = weight;
	array[pos++] = corrective;
	array[pos++] = trunc;
}

double EMTruncSqdist::Cluster::density_at(double phase) const {
	if (abs(phase - atan2(yy, xx)) < theta) {
		return trunc;
	} else {
		double dx = xx - cos(phase);
		double dy = yy - sin(phase);
		return corrective / (dx*dx + dy*dy);
	}
}

void EMTruncSqdist::Cluster::precalc() {
	double x = sqrt(xx*xx + yy*yy);
	double x2 = x*x;
	trunc = 1/(x2 + 1 - 2*x*cos(theta));
	corrective = 1/(2*(theta*trunc + (M_PI + 2*atan((x+1)*tan(theta/2)/(x-1)))/(1-x2)));
	trunc *= corrective;
}

void EMTruncSqdist::Cluster::estimate(double m1r, double m1i, double m2r, double m2i) {
	if (m1r == 0 && m2r == 0) {
		xx = 0;
		yy = 0;
		theta = 0;
		precalc();
		return;
	}

	int status;
	int tries = 0;
	double mu = atan2(m1i, m1r);
	const size_t n = 2;
	RightSide rs(m1r, m1i, m2r, m2i);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_multiroot_fsolver* s;
	const gsl_multiroot_fsolver_type* T;
	gsl_multiroot_function f = {&Cluster::estimation_f, n, &rs};

	do {
		int iter = 0;

		for (int i = 0; i < n; ++i) {
			gsl_vector_set(x, i, float(rand()) / RAND_MAX);
		}

		T = gsl_multiroot_fsolver_hybrids;
		s = gsl_multiroot_fsolver_alloc (T, n);
		gsl_multiroot_fsolver_set (s, &f, x);

		do
		{
			iter++;
			//std::cout << "x=(" << gsl_vector_get(s->x, 0) << ", " << gsl_vector_get(s->x, 1) << ") f=(" << gsl_vector_get(s->f, 0) << ", " << gsl_vector_get(s->f, 1) << ")" << std::endl;
			status = gsl_multiroot_fsolver_iterate(s);

			if (status) {
				break;
			}

			status = gsl_multiroot_test_residual (s->f, 1E-7);
		}	while (status == GSL_CONTINUE && iter < 1000);

		//std::cout << "status = " << gsl_strerror (status) << std::endl;
		tries++;
	} while (status != GSL_SUCCESS && tries < 1000);

	double xr = tanh(gsl_vector_get(s->x, 0));
	xr *= xr;
	//double xr = gsl_vector_get(s->x, 0);
	xx = xr*cos(mu);
	yy = xr*sin(mu);
	theta = tanh(gsl_vector_get(s->x, 1));
	theta *= M_PI*theta;
	//theta = gsl_vector_get(s->x, 1);

	if (theta < 0.05) theta = 0.05;

	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);

	precalc();
}

int EMTruncSqdist::Cluster::estimation_f(const gsl_vector* X, void* params, gsl_vector* f) {
	double x = tanh(gsl_vector_get(X, 0));
	x *= x;
	double t = tanh(gsl_vector_get(X, 1));
	t *= M_PI*t;
	/*double x = gsl_vector_get(X, 0);
	double t = gsl_vector_get(X, 1);*/
	RightSide* rs = (RightSide*) params;

	double x2 = x*x;
	double A = atan((x+1)*tan(t/2)/(x-1));
	double B = x2 + 1 - 2*x*cos(t);
	double C = t/A + (M_PI + 2*B)/(1 - x2);

	double m1 = (sin(t)/B + (M_PI*x + (x+1/x)*A)/(1 - x2) + t/(2*x))/C;
	double m2 = (sin(2*t)/(2*B) + (M_PI*x2 + (x2+1/x2)*A)/(1 - x2) + (t*x2 + 2*x*sin(t)+t)/(2*x2))/C;

	//std::cout << "(" << x << ", " << t << ") " << A << " " << B << " " << C << std::endl;

	gsl_vector_set(f, 0, m1 - rs->m1);
	gsl_vector_set(f, 1, m2 - rs->m2);

	return GSL_SUCCESS;
}

EMTruncSqdist::Cluster::RightSide::RightSide(double m1r, double m1i, double m2r, double m2i) :
	m1(sqrt(m1r*m1r + m1i*m1i)),
	m2(sqrt(m2r*m2r + m2i*m2i))
{

}

void EMTruncSqdist::save(FILE* file, bool lossy)
{
	int size = clusters.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		clusters[i].save(file, lossy);
	}
}

void EMTruncSqdist::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].load(file);
	}
}

void EMTruncSqdist::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = clusters.size();
	for (int i = 0; i < clusters.size(); ++i) {
		clusters[i].exportToArray(array, maxLen, pos);
	}
}

void EMTruncSqdist::importFromArray(double* array, int len, int& pos)
{
	int size = array[pos++];
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].importFromArray(array, len, pos);
	}
}

void EMTruncSqdist::print() {
	std::cout << "[";
	double sum_w = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		if (i) {
			std::cout << ", ";
		}
		clusters[i].print();
		sum_w += clusters[i].weight;
	}
	std::cout << "] (w = " << sum_w << ")" << std::endl;
}

double EMTruncSqdist::get_density_at(uint32_t time) {
	double phase = time_to_phase(time);
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(phase);
	}
	return result;
}
