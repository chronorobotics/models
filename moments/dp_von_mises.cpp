#include <gsl/gsl_sf_bessel.h>
#include <ctime>
#include "../CMoments.h"

#include "dp_von_mises.h"
#include "me_circular.h"

DPVonMises::DPVonMises(CMoments* parent_, int count_) :
	DensityParams(parent_, count_),
	kappa(),
	mu(),
	weight(),
	estimator()
{
	kappa.resize(count);
	mu.resize(count);
	weight.resize(count);
}

DPVonMises::~DPVonMises() {

}

MomentEstimator* DPVonMises::get_moment_estimator() {
	if (!estimator.get()) {
		estimator = std::shared_ptr<MECircular>(new MECircular(ceil(float(count)*3/2)));
	}
	return estimator.get();
}

void DPVonMises::reset(int count_) {
	count = count_;
	kappa.resize(count);
	mu.resize(count);
	weight.resize(count);
}

void DPVonMises::calculate()
{
	EquationParams ep(estimator->get_moments(), count, estimator->get_moment_count());
	ep.min_kappa = 1;
	const size_t n = ep.right_side.size();
	//int status;
	//int tries = 0;
	srandom(time(0));
	std::vector<double> buf1;
	std::vector<double> buf2;

	buf1.reserve(n);
	buf2.resize(n);

	for (int i = n-1; i >= 0; --i) {
		buf1.push_back(float(rand()) / RAND_MAX * 6);
	}

	std::vector<double>* result = fixed_point_iteration(&buf1, &buf2, ep);

	for (int i = 0; i < count; ++i) {
		kappa[i]  = lnhyp((*result)[3*i], ep.min_kappa);
		mu[i]     = (*result)[3*i + 1];
		weight[i] = hyp((*result)[3*i + 2]);
	}
}

double DPVonMises::density_at(uint32_t time) {
	double phase = MECircular::time_to_phase(time);
	double result = 0;
	for (int i = 0; i < count; ++i) {
		result += weight[i] * exp(kappa[i] * cos(phase - mu[i])) / (2 * M_PI * gsl_sf_bessel_I0(kappa[i]));
	}
	return result;
}

double DPVonMises::lnhyp(double x, double min_kappa) {
	//return log(1 + x/2 + sqrt(x*x + 4)/2);
	return log(x*x + 1) + min_kappa;
}

double DPVonMises::hyp(double x) {
	return sqrt(x*x + 1)-1;
}

double DPVonMises::dlnhyp(double x) {
	return 2*x / (x*x + 1);
}

double DPVonMises::dhyp(double x) {
	return x / sqrt(x*x + 1);
}

void DPVonMises::print() {
	std::cout << "[";
	for (int i = 0; i < count; ++i) {
		std::cout << "(" << kappa[i] << ", " << mu[i] << ", " << weight[i] << "), ";
	}
	std::cout << std::endl;
}

int DPVonMises::moment_f(const gsl_vector* x, void* params, gsl_vector* f) {
	EquationParams* ep = (EquationParams*) params;

	for (int i = ep->moment_count - 1; i >= 0; --i) {
		double y_re = 0;
		double y_im = 0;

		for (int j = ep->cluster_count - 1; j >= 0; --j) {
			double x_kappa  = lnhyp(gsl_vector_get(x, 3*j), ep->min_kappa);
			double x_mu     = gsl_vector_get(x, 3*j + 1);
			double x_weight = hyp(gsl_vector_get(x, 3*j + 2));

			if (isnan(x_kappa) || isnan(x_mu) || isnan(x_weight)) {
				return GSL_EDOM;
			}

			double foo;

			if (x_kappa > 500) {
				foo = 1;
			} else {
				foo = gsl_sf_bessel_In(i+1, x_kappa) / gsl_sf_bessel_I0(x_kappa);
			}

			y_re += x_weight * foo * cos((i+1)*x_mu);
			y_im += x_weight * foo * sin((i+1)*x_mu);
		}

		gsl_vector_set (f, 2*i    , y_re - ep->right_side[2*i]);
		if (2*i + 1 < ep->cluster_count*3) {
			gsl_vector_set (f, 2*i + 1, y_im - ep->right_side[2*i + 1]);
		}

	}

	return GSL_SUCCESS;
}

int DPVonMises::save(FILE* file, bool lossy)
{
	Distribution type = VON_MISES;
	fwrite(&type, sizeof(Distribution), 1, file);

	fwrite(&count, sizeof(int), 1, file);
	fwrite(&kappa[0], sizeof(double), count, file);
	fwrite(&mu[0], sizeof(double), count, file);
	fwrite(&weight[0], sizeof(double), count, file);
	return 0;
}

int DPVonMises::load(FILE* file)
{
	int cnt;
	fread(&cnt, sizeof(int), 1, file);
	reset(cnt);
	fread(&kappa[0], sizeof(double), cnt, file);
	fread(&mu[0], sizeof(double), cnt, file);
	fread(&weight[0], sizeof(double), cnt, file);
	return 0;
}


void DPVonMises::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = (double) VON_MISES;
	array[pos++] = count;
	for (int i = 0; i <= count; ++i) {
		array[pos++] = kappa[i];
		array[pos++] = mu[i];
		array[pos++] = weight[i];
	}
}

void DPVonMises::importFromArray(double* array, int len, int& pos)
{
	int cnt;
	cnt = array[pos++];
	reset(cnt);
	for (int i = 0; i <= count; ++i) {
		kappa[i] = array[pos++];
		mu[i] = array[pos++];
		weight[i] = array[pos++];
	}
}

DPVonMises::EquationParams::EquationParams(std::vector<double> right_side_, int cluster_count_, int moment_count_) :
	right_side(right_side_),
	cluster_count(cluster_count_),
	moment_count(moment_count_),
	min_kappa()
{

}

double DPVonMises::InI0(int n, double kappa) const {
	return gsl_sf_bessel_In(n, kappa) / gsl_sf_bessel_I0(kappa);
}

double DPVonMises::dInI0(int n, double kappa) const {
	double i0 = gsl_sf_bessel_I0(kappa);
	return (i0*(gsl_sf_bessel_In(n-1, kappa) + gsl_sf_bessel_In(n+1, kappa)) + 2*gsl_sf_bessel_In(n, kappa)*gsl_sf_bessel_I1(kappa)) / (2*i0*i0);
}

double DPVonMises::segsgn(double x) const {
	double result = sqrt(1 + x*x);
	if (x < 0) {
		return -result;
	} else {
		return result;
	}
}

double DPVonMises::fpi_d(int eq, std::vector<double>* now, const EquationParams& ep) const {
	int i = eq / 2;
	int j = eq / 3;

	double x_kappa  = lnhyp((*now)[3*j], ep.min_kappa);
	double x_mu     = (*now)[3*j + 1];
	double x_weight = (*now)[3*j + 2];

	switch (eq % 3) {
		case 0: {
			double foo = dInI0(i+1, x_kappa) * hyp(x_weight) * dlnhyp((*now)[3*j]);
			if (i % 2) {
				return sin((i+1)*x_mu) * foo;
			} else {
				return cos((i+1)*x_mu) * foo;
			}
		} break;
		case 1: {
			double foo = InI0(i+1, x_kappa) * hyp(x_weight) * (i+1);
			if (i % 2) {
				return cos((i+1)*x_mu) * foo;
			} else {
				return -sin((i+1)*x_mu) * foo;
			}
		} break;
		case 2: {
			double foo = InI0(i+1, x_kappa) * dhyp(x_weight);
			if (i % 2) {
				return sin((i+1)*x_mu) * foo;
			} else {
				return cos((i+1)*x_mu) * foo;
			}
		} break;
		default:
			return 0;
			break;
	}
}

std::vector<double>* DPVonMises::fixed_point_iteration(std::vector<double>* now, std::vector<double>* next, const EquationParams& ep) const {
	int iter = 0;
	float dist = 10000;
	while (dist > 0.1) {

		std::cout << "iter " << iter << std::endl;
		std::cout << std::endl << "    value:";
		for (int i = 0; i < now->size(); ++i) {
			std::cout << " " << (*now)[i];
		}
		std::cout << std::endl;

		for (int i = ep.moment_count - 1; i >= 0; --i) {
			double y_re = 0;
			double y_im = 0;

			for (int j = ep.cluster_count - 1; j >= 0; --j) {
				double x_kappa  = lnhyp((*now)[3*j], ep.min_kappa);
				double x_mu     = (*now)[3*j + 1];
				double x_weight = hyp((*now)[3*j + 2]);

				double foo = InI0(i+1, x_kappa);
				double sx = sin((i+1)*x_mu);
				double cx = cos((i+1)*x_mu);

				y_re += x_weight * foo * cx;
				y_im += x_weight * foo * sx;
			}

			int var1 = 2*i;
			int var2 = 2*i + 1;

			y_re = (y_re - ep.right_side[var1]) / (i+1);
			y_im = (y_im - ep.right_side[var2]) / (i+1);


			(*next)[var1] = (*now)[var1] - y_re/segsgn(fpi_d(var1, now, ep));
			if (var2 < ep.cluster_count*3) {
				(*next)[var2] = (*now)[var2] - y_im/segsgn(fpi_d(var2, now, ep));
			}

		}

		dist = 0;
		for (int i = ep.cluster_count*3-1; i >= 0; --i) {
			dist += ((*now)[i]-(*next)[i]) * ((*now)[i]-(*next)[i]);
		}
		dist = sqrt(dist);

		std::swap(now, next);
		iter++;
	}

	return next;
}
