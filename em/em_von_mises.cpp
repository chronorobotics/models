#include <iostream>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include "em_von_mises.h"

EMVonMises::EMVonMises() :
	EMCircular(),
	clusters(),
	timestamps_weight(0)
{

}

EMVonMises::EMVonMises(int cluster_count_) :
	EMCircular(cluster_count_),
	clusters(),
	timestamps_weight(0)
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

	for (unsigned int i = 0; i < timestamps.size(); ++i) {
		double s = 0;
		for (unsigned int j = 0; j < clusters.size(); ++j) {
			double a = clusters[j].density_at(timestamps[i].phase) * clusters[j].weight;
			timestamps[i].alpha[j] = a;
			s += a;
		}

		for (unsigned int j = 0; j < clusters.size(); ++j) {
			timestamps[i].alpha[j] /= s;
		}
	}
}

double EMVonMises::maximisation(bool keep_kappa) {
	double shift = 0;
	std::cerr << "Performing maximisation ..." << std::endl;

	for (unsigned int i = 0; i < clusters.size(); ++i) {
		double last_kappa = clusters[i].kappa;
		double last_mu = clusters[i].mu;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (unsigned int j = 0; j < timestamps.size(); ++j) {
			s += timestamps[j].alpha[i];
		}
		clusters[i].weight = s / timestamps.size();

		double mean_re = 0;
		double mean_im = 0;
		for (unsigned int j = 0; j < timestamps.size(); ++j) {
			mean_re += cos(timestamps[j].phase) * timestamps[j].alpha[i];
			mean_im += sin(timestamps[j].phase) * timestamps[j].alpha[i];
		}
		mean_re /= s;
		mean_im /= s;
		clusters[i].estimate_from_mean(mean_re, mean_im, keep_kappa);

		double delta_kappa = last_kappa - clusters[i].kappa;
		double delta_mu = last_mu - clusters[i].mu;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_kappa*delta_kappa + delta_mu*delta_mu + delta_weight*delta_weight;
	}

	return sqrt(shift);
}

void EMVonMises::train() {
	double shift = 555;
	double dl;
	do {
		double l = get_loglikelihood();
		expectation();
		shift = maximisation(false);
		dl = get_loglikelihood() - l;
		std::cout << "dl = " << dl << std::endl;
	} while (isnan(shift) || dl > 1.0);
}

void EMVonMises::add_time(uint32_t time, double value) {
	timestamps.push_back(Timestamp(time, cluster_count, value));
	timestamps_weight += value;
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

void EMVonMises::Cluster::estimate_from_mean(double re, double im, bool keep_kappa) {
	if (re == 0 && im == 0) {
		mu = 0;
		kappa = 0;
		return;
	}

	mu = atan2(im, re);
	if (keep_kappa) {
		return;
	}

	double aa = sqrt(re*re + im*im);

	int status;
	int iter = 0;
	int max_iter = 1000;
	double x0;
	double x = 5.0;

	if (mean_f(300, &aa) < 0) {
		kappa = 300;
		return;
	}

	const gsl_root_fsolver_type* T = gsl_root_fsolver_bisection;
	gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);

	gsl_function FDF;
	FDF.function = &mean_f;
	//FDF.df = &mean_df;
	//FDF.fdf = &mean_fdf;
	FDF.params = &aa;

	gsl_root_fsolver_set (s, &FDF, 0, 300);

	do {
		iter++;
		status = gsl_root_fsolver_iterate (s);
		x0 = x;
		x = gsl_root_fsolver_root (s);
		status = gsl_root_test_delta (x, x0, 0, 1e-30);

	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
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

double EMVonMises::get_density_at_d(double phase) const {
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(phase);
	}
	return result;
}
