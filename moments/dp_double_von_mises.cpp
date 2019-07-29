#include <gsl/gsl_sf_bessel.h>
#include <ctime>
#include "../CMoments.h"

#include "dp_double_von_mises.h"

DPDoubleVonMises::DPDoubleVonMises(CMoments* parent_, int count_) :
	DensityParams(parent_, count_),
	kappa(),
	mu(),
	lambda(),
	nu(),
	weight(),
	estimator()
{
	kappa.resize(count);
	mu.resize(count);
	lambda.resize(count);
	nu.resize(count);
	weight.resize(count);
}

DPDoubleVonMises::~DPDoubleVonMises() {

}

MomentEstimator* DPDoubleVonMises::get_moment_estimator() {
	if (!estimator.get()) {
		estimator = std::shared_ptr<MEBicircular>(new MEBicircular(ceil(float(count)*5/2)));
	}
	return estimator.get();
}

void DPDoubleVonMises::reset(int count_) {
	count = count_;
	kappa.resize(count);
	mu.resize(count);
	lambda.resize(count);
	nu.resize(count);
	weight.resize(count);
}

void DPDoubleVonMises::calculate()
{
	EquationParams ep(estimator->get_moments(), estimator->get_moment_indices(), count, estimator->get_moment_count());
	int status;
	int tries = 0;
	srandom(time(0));

	do {
		tries++;
		ep.min_kappa = tan(float(rand()) / RAND_MAX * M_PI_2);
		ep.min_lambda = tan(float(rand()) / RAND_MAX * M_PI_2);
		int iter = 0;
		const size_t n = ep.right_side.size();
		gsl_multiroot_function f = {&DPDoubleVonMises::moment_f, n, &ep};

		gsl_vector* x = gsl_vector_alloc(n);
		for (int i = 0; i < n; ++i) {
			gsl_vector_set(x, i, float(rand()) / RAND_MAX * 6);
		}

		gsl_multiroot_fsolver* s;
		const gsl_multiroot_fsolver_type* T;
		T = gsl_multiroot_fsolver_hybrids;
		s = gsl_multiroot_fsolver_alloc (T, n);
		gsl_multiroot_fsolver_set (s, &f, x);

		do
		{
			iter++;
			status = gsl_multiroot_fsolver_iterate(s);

			std::cout << "iter " << iter << std::endl;
			std::cout << "    residuum:";
			for (int i = 0; i < n; ++i) {
				std::cout << " " << gsl_vector_get(s->f, i);
			}
			std::cout << std::endl << "    value:";
			for (int i = 0; i < n; ++i) {
				std::cout << " " << gsl_vector_get(s->x, i);
			}
			std::cout << std::endl;

			if (status) {
				break;
			}

			status = gsl_multiroot_test_residual (s->f, 0.1);
		}	while (status == GSL_CONTINUE && iter < 1000);

		std::cout << /*"\33[2K\" <<*/ "rstatus = " << gsl_strerror (status) << ", tries = " << tries << std::endl;
		for (int i = 0; i < count; ++i) {
			kappa[i]  = lnhyp(gsl_vector_get(s->x, 5*i), ep.min_kappa);
			mu[i]     = gsl_vector_get(s->x, 5*i + 1);
			lambda[i]  = lnhyp(gsl_vector_get(s->x, 5*i + 2), ep.min_lambda);
			nu[i]     = gsl_vector_get(s->x, 5*i + 3);
			weight[i] = hyp(gsl_vector_get(s->x, 5*i + 4));
		}

		gsl_multiroot_fsolver_free(s);
		gsl_vector_free(x);
	} while (status != GSL_SUCCESS && tries < 1);

	std::cout << std::endl;

	if (status != GSL_SUCCESS) {
		std::cout << "ERROR: Solution not found!" << std::endl;
	}
}

double DPDoubleVonMises::density_at(uint32_t time) {
	MEBicircular::Phase phase(time);
	double result = 0;
	for (int i = 0; i < count; ++i) {
		result += weight[i] * exp(kappa[i] * cos(phase.f1 - mu[i])) / (2 * M_PI * gsl_sf_bessel_I0(kappa[i]))
							* exp(lambda[i] * cos(phase.f2 - nu[i])) / (2 * M_PI * gsl_sf_bessel_I0(lambda[i]));
	}
	return result;
}

double DPDoubleVonMises::lnhyp(double x, double min_kappa) {
	return log(x*x + 1) + min_kappa;
}

double DPDoubleVonMises::hyp(double x) {
	return sqrt(x*x + 1)-1;
}

void DPDoubleVonMises::print() {
	std::cout << "[";
	for (int i = 0; i < count; ++i) {
		std::cout << "(" << kappa[i] << ", " << mu[i] << ", " << lambda[i] << ", " << nu[i] << ", " << weight[i] << "), ";
	}
	std::cout << std::endl;
}

int DPDoubleVonMises::moment_f(const gsl_vector* x, void* params, gsl_vector* f) {
	EquationParams* ep = (EquationParams*) params;

	for (int i = ep->moment_count - 1; i >= 0; --i) {
		double y_re = 0;
		double y_im = 0;

		for (int j = ep->cluster_count - 1; j >= 0; --j) {
			double x_kappa  = lnhyp(gsl_vector_get(x, 5*j), ep->min_kappa);
			double x_mu     = gsl_vector_get(x, 5*j + 1);
			double x_lambda  = lnhyp(gsl_vector_get(x, 5*j + 2), ep->min_lambda);
			double x_nu     = gsl_vector_get(x, 5*j + 3);
			double x_weight = hyp(gsl_vector_get(x, 5*j + 4));

			if (isnan(x_kappa) || isnan(x_mu) || isnan(x_lambda) || isnan(x_nu) || isnan(x_weight)) {
				return GSL_EDOM;
			}

			double foo;

			if (x_kappa > 500) {
				foo = 1;
			} else {
				foo = gsl_sf_bessel_In(ep->indices[i].i1, x_kappa) * gsl_sf_bessel_In(ep->indices[i].i2, x_lambda)
						/ (gsl_sf_bessel_I0(x_kappa) * gsl_sf_bessel_I0(x_lambda));
			}

			y_re += x_weight * foo * cos(ep->indices[i].i1*x_mu) * cos(ep->indices[i].i2*x_nu);
			y_im += x_weight * foo * sin(ep->indices[i].i1*x_mu) * sin(ep->indices[i].i2*x_nu);
		}

		gsl_vector_set (f, 2*i    , y_re - ep->right_side[2*i]);
		if (2*i + 1 < ep->cluster_count*5) {
			gsl_vector_set (f, 2*i + 1, y_im - ep->right_side[2*i + 1]);
		}

	}

	return GSL_SUCCESS;
}

int DPDoubleVonMises::save(FILE* file, bool lossy)
{
	Distribution type = VON_MISES;
	fwrite(&type, sizeof(Distribution), 1, file);

	fwrite(&count, sizeof(int), 1, file);
	fwrite(&kappa[0], sizeof(double), count, file);
	fwrite(&mu[0], sizeof(double), count, file);
	fwrite(&lambda[0], sizeof(double), count, file);
	fwrite(&nu[0], sizeof(double), count, file);
	fwrite(&weight[0], sizeof(double), count, file);
	return 0;
}

int DPDoubleVonMises::load(FILE* file)
{
	int cnt;
	fread(&cnt, sizeof(int), 1, file);
	reset(cnt);
	fread(&kappa[0], sizeof(double), cnt, file);
	fread(&mu[0], sizeof(double), cnt, file);
	fread(&lambda[0], sizeof(double), cnt, file);
	fread(&nu[0], sizeof(double), cnt, file);
	fread(&weight[0], sizeof(double), cnt, file);
	return 0;
}


void DPDoubleVonMises::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = (double) VON_MISES;
	array[pos++] = count;
	for (int i = 0; i <= count; ++i) {
		array[pos++] = kappa[i];
		array[pos++] = mu[i];
		array[pos++] = lambda[i];
		array[pos++] = nu[i];
		array[pos++] = weight[i];
	}
}

void DPDoubleVonMises::importFromArray(double* array, int len, int& pos)
{
	int cnt;
	cnt = array[pos++];
	reset(cnt);
	for (int i = 0; i <= count; ++i) {
		kappa[i] = array[pos++];
		mu[i] = array[pos++];
		lambda[i] = array[pos++];
		nu[i] = array[pos++];
		weight[i] = array[pos++];
	}
}

DPDoubleVonMises::EquationParams::EquationParams(std::vector<double> right_side_, const std::vector<MEBicircular::Index>& indices_,
																								 int cluster_count_, int moment_count_) :
	right_side(right_side_),
	indices(indices_),
	cluster_count(cluster_count_),
	moment_count(moment_count_),
	min_kappa(),
	min_lambda()
{

}
