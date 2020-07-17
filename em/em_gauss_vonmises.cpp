#include <math.h>
#include <iostream>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>

#include "em_gauss_vonmises.h"

EMGaussVonMises::EMGaussVonMises() :
	cluster_count(),
	total_weight(0),
	clusters(),
	points()
{

}

EMGaussVonMises::EMGaussVonMises(int cluster_count_) :
	cluster_count(cluster_count_),
	total_weight(0),
	clusters(),
	points()
{
	double s = 0;
	for (int i = 0; i < cluster_count; ++i) {
		clusters.push_back(Cluster());
		s += clusters[i].weight;
	}

	for (int i = 0; i < cluster_count; ++i) {
		clusters[i].weight /= s;
	}
}

double EMGaussVonMises::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 86400.0f) / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

void EMGaussVonMises::expectation() {
	std::cerr << "Performing expectation ..." << std::endl;

	for (int i = 0; i < points.size(); ++i) {
		double s = 0;
		for (int j = 0; j < clusters.size(); ++j) {
			double a = clusters[j].density_at(points[i].x, points[i].y, points[i].phase);
			points[i].alpha[j] = a;
			s += a;
		}

		for (int j = 0; j < clusters.size(); ++j) {
			points[i].alpha[j] /= s;
		}
	}
}

double EMGaussVonMises::maximisation() {
	double shift = 0;
	std::cerr << "Performing maximisation ..." << std::endl;
	for (int i = 0; i < clusters.size(); ++i) {
		double last_ex = clusters[i].ex;
		double last_ey = clusters[i].ey;
		double last_cxx = clusters[i].cxx;
		double last_cyy = clusters[i].cyy;
		double last_cxy = clusters[i].cxy;
		double last_tmu = clusters[i].tmu;
		double last_tkappa = clusters[i].tkappa;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < points.size(); ++j) {
			s += points[j].alpha[i];
		}
		clusters[i].weight = s / points.size();

		double mean_x = 0;
		double mean_y = 0;
		for (int j = 0; j < points.size(); ++j) {
			mean_x += points[j].x * points[j].alpha[i];
			mean_y += points[j].y * points[j].alpha[i];
		}
		mean_x /= s;
		mean_y /= s;
		clusters[i].ex = mean_x;
		clusters[i].ey = mean_y;

		double mean_re = 0;
		double mean_im = 0;
		for (unsigned int j = 0; j < points.size(); ++j) {
			mean_re += cos(points[j].phase) * points[j].alpha[i];
			mean_im += sin(points[j].phase) * points[j].alpha[i];
		}
		mean_re /= s;
		mean_im /= s;
		clusters[i].estimate_from_mean(mean_re, mean_im, false);

		double deviation_x = 0;
		double deviation_y = 0;
		double covariance = 0;
		for (int j = 0; j < points.size(); ++j) {
			double dx = points[j].x - mean_x;
			double dy = points[j].y - mean_y;
			deviation_x += dx*dx * points[j].alpha[i];
			deviation_y += dy*dy * points[j].alpha[i];
			covariance += dx*dy * points[j].alpha[i];
		}
		deviation_x /= s;
		deviation_y /= s;
		covariance /= s;
		double foo = sqrt(deviation_x*deviation_y - covariance*covariance);
		if (foo <= 0.4) {
			foo = 0.4 / foo;
			deviation_x *= foo;
			deviation_y *= foo;
			covariance *= foo;
		}

		clusters[i].det = deviation_x*deviation_y - covariance*covariance;
		clusters[i].cxx = deviation_y / clusters[i].det;
		clusters[i].cyy = deviation_x / clusters[i].det;
		clusters[i].cxy = -2*covariance / clusters[i].det;
		clusters[i].det = sqrt(4*M_PI*M_PI*clusters[i].det);

		double delta_ex = last_ex - clusters[i].ex;
		double delta_ey = last_ey - clusters[i].ey;
		double delta_cxx = last_cxx - clusters[i].cxx;
		double delta_cyy = last_cyy - clusters[i].cyy;
		double delta_cxy = last_cxy - clusters[i].cxy;
		double delta_tmu = last_tmu - clusters[i].tmu;
		double delta_tkappa = last_tkappa - clusters[i].tkappa;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_ex*delta_ex + delta_ey*delta_ey + delta_cxx*delta_cxx + delta_cyy*delta_cyy +
						 delta_cxy*delta_cxy + delta_tmu*delta_tmu + delta_tkappa*delta_tkappa +
						 delta_weight*delta_weight;
	}
	return sqrt(shift);
}

void EMGaussVonMises::train() {
	double shift = 555;
	do {
		expectation();
		shift = maximisation();
		std::cerr << "shift = " << shift << std::endl;
	} while (shift > 1E-3);
}

void EMGaussVonMises::add_value(double x, double y, uint32_t time, double weight) {
	points.push_back(Point(x, y, time, cluster_count, weight));
	total_weight += weight;
}

EMGaussVonMises::Cluster::Cluster() :
	ex(2*float(rand()) / RAND_MAX - 1),
	ey(2*float(rand()) / RAND_MAX - 1),
	cxx(),
	cyy(),
	cxy(),
	det(),
	tmu(float(rand()) / RAND_MAX),
	tkappa(3 + float(rand())/RAND_MAX),
	weight(float(rand()) / RAND_MAX)
{
	double a = float(rand()) / RAND_MAX + 1;
	double b = float(rand()) / RAND_MAX + 1;
	double c = float(rand()) / RAND_MAX * sqrt(a*b);
	det = a*b - c*c;
	cxx = b / det;
	cyy = a / det;
	cxy = -2*c / det;
	det = sqrt(4*M_PI*M_PI*det);
}

EMGaussVonMises::Cluster::Cluster(double ex_, double ey_, double cxx_, double cyy_,
																double cxy_, double det_, double tmu_, double tkappa_,
																double weight_) :
	ex(ex_),
	ey(ey_),
	cxx(cxx_),
	cyy(cyy_),
	cxy(cxy_),
	det(det_),
	tmu(tmu_),
	tkappa(tkappa_),
	weight(weight_)
{

}

void EMGaussVonMises::Cluster::print() {
	double det_ = det*det/4/M_PI/M_PI;
	double DX = cyy*det_;
	double DY = cxx*det_;
	double XY = -cxy*det_/2;
	std::cout << "(" << ex << ", " << ey << ", " << DX << ", " << DY << ", " << XY << ", " << tmu << ", " << tkappa << ", " << weight << ")";
}

void EMGaussVonMises::Cluster::save(FILE* file, bool lossy) {
	fwrite(&ex, sizeof(double), 1, file);
	fwrite(&ey, sizeof(double), 1, file);
	fwrite(&cxx, sizeof(double), 1, file);
	fwrite(&cxy, sizeof(double), 1, file);
	fwrite(&cyy, sizeof(double), 1, file);
	fwrite(&det, sizeof(double), 1, file);
	fwrite(&tmu, sizeof(double), 1, file);
	fwrite(&tkappa, sizeof(double), 1, file);
	fwrite(&weight, sizeof(double), 1, file);
}

void EMGaussVonMises::Cluster::load(FILE* file) {
	fread(&ex, sizeof(double), 1, file);
	fread(&ey, sizeof(double), 1, file);
	fread(&cxx, sizeof(double), 1, file);
	fread(&cxy, sizeof(double), 1, file);
	fread(&cyy, sizeof(double), 1, file);
	fread(&det, sizeof(double), 1, file);
	fread(&tmu, sizeof(double), 1, file);
	fread(&tkappa, sizeof(double), 1, file);
	fread(&weight, sizeof(double), 1, file);
}

void EMGaussVonMises::Cluster::importFromArray(double* array, int len, int& pos) {
	ex = array[pos++];
	ey = array[pos++];
	cxx = array[pos++];
	cyy = array[pos++];
	cxy = array[pos++];
	det = array[pos++];
	tmu = array[pos++];
	tkappa = array[pos++];
	weight = array[pos++];
}

void EMGaussVonMises::Cluster::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = ex;
	array[pos++] = ey;
	array[pos++] = cxx;
	array[pos++] = cyy;
	array[pos++] = cxy;
	array[pos++] = det;
	array[pos++] = tmu;
	array[pos++] = tkappa;
	array[pos++] = weight;
}

double EMGaussVonMises::Cluster::density_at(double x, double y, double phase) const {
	double VonMises = exp(tkappa*cos(phase - tmu)) / (2*M_PI*gsl_sf_bessel_I0(tkappa));
	double dx = ex - x;
	double dy = ey - y;
	double gauss = exp(-(cxx*dx*dx + cyy*dy*dy + cxy*dx*dy)/2) / det;
	return VonMises * gauss;
}

void EMGaussVonMises::save(FILE* file, bool lossy)
{
	int size = clusters.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		clusters[i].save(file, lossy);
	}
}

void EMGaussVonMises::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].load(file);
	}
}

void EMGaussVonMises::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = clusters.size();
	for (int i = 0; i < clusters.size(); ++i) {
		clusters[i].exportToArray(array, maxLen, pos);
	}
}

void EMGaussVonMises::importFromArray(double* array, int len, int& pos)
{
	int size = array[pos++];
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].importFromArray(array, len, pos);
	}
}

void EMGaussVonMises::print() {
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

double EMGaussVonMises::get_density_at(double x, double y, uint32_t time) const {
	double phase = time_to_phase(time);
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(x, y, phase);
	}
	return result;
}

EMGaussVonMises::Point::Point(double x_, double y_, uint32_t time, int cluster_count, double weight_) :
	x(x_),
	y(y_),
	phase(time_to_phase(time)),
	alpha(),
	weight(weight_)
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

void EMGaussVonMises::Cluster::estimate_from_mean(double re, double im, bool keep_kappa) {
	if (re == 0 && im == 0) {
		tmu = 0;
		tkappa = 0;
		return;
	}

	tmu = atan2(im, re);
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
		tkappa = 300;
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
	tkappa = x;
}

double EMGaussVonMises::Cluster::mean_f (double x, void* params) {
	double* y = (double*) params;
	return gsl_sf_bessel_I1(x) / gsl_sf_bessel_I0(x) - *y;
}

double EMGaussVonMises::Cluster::mean_df (double x, void* params) {
	double i0 = gsl_sf_bessel_I0(x);
	double i1 = gsl_sf_bessel_I1(x);
	return (i0*(i0 + gsl_sf_bessel_In(2, x)) + 2*i1*i1) / (i0*i0);
}

void EMGaussVonMises::Cluster::mean_fdf (double x, void* params, double* y, double* dy) {
	double* y_ = (double*) params;
	double i0 = gsl_sf_bessel_I0(x);
	double i1 = gsl_sf_bessel_I1(x);

	*y = i1/i0 - *y_;
	*dy = (i0*(i0 + gsl_sf_bessel_In(2, x)) + 2*i1*i1) / (i0*i0);
}
