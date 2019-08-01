#include <iostream>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include "em_von_mises.h"

EMVonMises::EMVonMises(int cluster_count_) :
	cluster_count(cluster_count_),
	clusters(),
	timestamps()
{
	double s = 0;
	for (int i = 0; i < cluster_count; ++i) {
		clusters.push_back(Cluster());
		s += clusters[i].weight;
		clusters[i].kappa = 374 + 146*float(rand())/RAND_MAX;
		clusters[i].mu = float(rand())/RAND_MAX * 2 * M_PI;
	}

	for (int i = 0; i < cluster_count; ++i) {
		clusters[i].weight /= s;
	}
}

void EMVonMises::expectation() {
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

double EMVonMises::maximisation() {
	double shift;
	std::cerr << "Performing maximisation ..." << std::endl;

	for (int i = 0; i < clusters.size(); ++i) {
		double last_kappa = clusters[i].kappa;
		double last_mu = clusters[i].mu;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			s += timestamps[j].alpha[i];
		}
		clusters[i].weight = s / timestamps.size();

		double mean_re = 0;
		double mean_im = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			mean_re += cos(timestamps[j].phase) * timestamps[j].alpha[i];
			mean_im += sin(timestamps[j].phase) * timestamps[j].alpha[i];
		}
		mean_re /= s;
		mean_im /= s;
		clusters[i].estimate_from_mean(mean_re, mean_im);

		double delta_kappa = last_kappa - clusters[i].kappa;
		double delta_mu = last_mu - clusters[i].mu;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_kappa*delta_kappa + delta_mu*delta_mu + delta_weight*delta_weight;
	}

	return sqrt(shift);
}

double EMVonMises::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 86400.0f) / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

void EMVonMises::train() {
	double shift;
	do {
		expectation();
		shift = maximisation();
		std::cout << "shift = " << shift << std::endl;
	} while (isnan(shift) || shift > 1E-7);
}

void EMVonMises::add_time(uint32_t time) {
	timestamps.push_back(Timestamp(time, cluster_count));
}

EMVonMises::Timestamp::Timestamp(uint32_t time, int cluster_count) :
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

EMVonMises::Cluster::Cluster() :
	kappa(3 + float(rand())/RAND_MAX),
	mu(float(rand()) / RAND_MAX),
	weight(float(rand()) / RAND_MAX)
{

}

EMVonMises::Cluster::Cluster(double kappa_, double mu_, double weight_) :
	kappa(kappa_),
	mu(mu_),
	weight()
{

}

void EMVonMises::Cluster::print() {
	std::cout << "(" << kappa << ", " << mu << ", " << weight << ")";
}

void EMVonMises::Cluster::save(FILE* file, bool lossy) {
	fwrite(&kappa, sizeof(double), 1, file);
	fwrite(&mu, sizeof(double), 1, file);
	fwrite(&weight, sizeof(double), 1, file);
}

void EMVonMises::Cluster::load(FILE* file) {
	fread(&kappa, sizeof(double), 1, file);
	fread(&mu, sizeof(double), 1, file);
	fread(&weight, sizeof(double), 1, file);
}

void EMVonMises::Cluster::importFromArray(double* array, int len, int& pos) {
	kappa = array[pos++];
	mu = array[pos++];
	weight = array[pos++];
}

void EMVonMises::Cluster::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = kappa;
	array[pos++] = mu;
	array[pos++] = weight;
}

double EMVonMises::Cluster::density_at(double phase) const {
	return exp(kappa*cos(phase - mu)) / (2*M_PI*gsl_sf_bessel_I0(kappa));
}

void EMVonMises::Cluster::estimate_from_mean(double re, double im) {
	if (re == 0 && im == 0) {
		mu = 0;
		kappa = 0;
		return;
	}

	mu = atan2(im, re);
	//return;
	double aa = sqrt(re*re + im*im);

	int status;
	int iter = 0;
	int max_iter = 1000;
	double x0;
	double x = 5.0;

	const gsl_root_fdfsolver_type* T = gsl_root_fdfsolver_newton;
	gsl_root_fdfsolver* s = gsl_root_fdfsolver_alloc(T);

	gsl_function_fdf FDF;
	FDF.f = &mean_f;
	FDF.df = &mean_df;
	FDF.fdf = &mean_fdf;
	FDF.params = &aa;

	gsl_root_fdfsolver_set (s, &FDF, x);

	//printf ("using %s method\n", gsl_root_fdfsolver_name (s));
	//printf ("%-5s %10s %10s %10s\n", "iter", "root", "err", "err(est)");

	do {
		iter++;
		status = gsl_root_fdfsolver_iterate (s);
		x0 = x;
		x = gsl_root_fdfsolver_root (s);
		status = gsl_root_test_delta (x, x0, 0, 1e-3);

		//if (status == GSL_SUCCESS)
			//std::cout << "Converged, iter=" << iter << std::endl;

		//printf ("%5d %10.7f %10.7f\n", iter, x, x - x0);
	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fdfsolver_free (s);
	kappa = x;
}

double EMVonMises::Cluster::mean_f (double x, void* params) {
	double* y = (double*) params;
	return gsl_sf_bessel_I1(x) / gsl_sf_bessel_I0(x) - *y;
}

double EMVonMises::Cluster::mean_df (double x, void* params) {
	double i0 = gsl_sf_bessel_I0(x);
	double i1 = gsl_sf_bessel_I1(x);
	return (i0*(i0 + gsl_sf_bessel_In(2, x)) + 2*i1*i1) / (i0*i0);
}

void EMVonMises::Cluster::mean_fdf (double x, void* params, double* y, double* dy) {
	double* y_ = (double*) params;
	double i0 = gsl_sf_bessel_I0(x);
	double i1 = gsl_sf_bessel_I1(x);

	*y = i1/i0 - *y_;
	*dy = (i0*(i0 + gsl_sf_bessel_In(2, x)) + 2*i1*i1) / (i0*i0);
}

void EMVonMises::save(FILE* file, bool lossy)
{
	int size = clusters.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		clusters[i].save(file, lossy);
	}
}

void EMVonMises::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].load(file);
	}
}

void EMVonMises::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = clusters.size();
	for (int i = 0; i < clusters.size(); ++i) {
		clusters[i].exportToArray(array, maxLen, pos);
	}
}

void EMVonMises::importFromArray(double* array, int len, int& pos)
{
	int size = array[pos++];
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].importFromArray(array, len, pos);
	}
}

void EMVonMises::print() {
	std::cout << "[";
	for (int i = 0; i < clusters.size(); ++i) {
		if (i) {
			std::cout << ", ";
		}
		clusters[i].print();
	}
	std::cout << "]" << std::endl;
}

double EMVonMises::get_density_at(uint32_t time) {
	double phase = time_to_phase(time);
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(phase);
	}
	return result;
}
