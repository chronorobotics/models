#include <gsl/gsl_sf_bessel.h>
#include "../CMoments.h"

#include "dp_von_mises.h"
#include "right_side.h"

DPVonMises::DPVonMises(int count_) :
	DensityParams(count_),
	kappa(),
	mu(),
	weight()
{
	if (count_ < 0) {
		count = CMoments::moment_count*2/3;
	} else {
		count = count_;
	}
	kappa = new double[count];
	mu = new double[count];
	weight = new double[count];
}

DPVonMises::~DPVonMises() {
	delete[] kappa;
	delete[] mu;
	delete[] weight;
}

void DPVonMises::reset(int count_) {
	delete[] kappa;
	delete[] mu;
	delete[] weight;
	count = count_;
	kappa = new double[count];
	mu = new double[count];
	weight = new double[count];
}

void DPVonMises::calculate(RightSide& rs)
{
	int status;
	int iter = 0;

	const size_t n = CMoments::moment_count*2;
	std::cout << n << " " << count << std::endl;
	gsl_multiroot_function f = {&DPVonMises::moment_f, n, &rs};

	gsl_vector *x = gsl_vector_alloc(n);
	for (int i = 0; i < n; ++i) {
		gsl_vector_set(x, i, 1.012431 + i*0.10347);
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

		std::cout << "iter "<< iter <<" "<< gsl_vector_get(s->f, 0) << " " << gsl_vector_get(s->f, 1) << " " << gsl_vector_get(s->f, 2) << " "
							<< gsl_vector_get(s->f, 3) << " " << gsl_vector_get(s->f, 4) << " " << gsl_vector_get(s->f, 5) << std::endl;
		std::cout << gsl_vector_get(s->x, 0) << " " << gsl_vector_get(s->x, 1) << " " << gsl_vector_get(s->x, 2) << " "
							<< gsl_vector_get(s->x, 3) << " " << gsl_vector_get(s->x, 4) << " " << gsl_vector_get(s->x, 5) << std::endl;

		if (status) {
			break;
		}

		status = gsl_multiroot_test_residual (s->f, 0.07);
	}	while (status == GSL_CONTINUE && iter < 1000);

	printf ("status = %s\n", gsl_strerror (status));
	for (int i = 0; i < count; ++i) {
		kappa[i]  = lnhyp(gsl_vector_get(s->x, 3*i));
		mu[i]     = gsl_vector_get(s->x, 3*i + 1) + M_PI_2;
		weight[i] = hyp(gsl_vector_get(s->x, 3*i + 2));
	}

	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
}

double DPVonMises::density_at(double phase) {
	double result = 0;
	for (int i = 0; i < count; ++i) {
		result += weight[i] * exp(kappa[i] * cos(phase - mu[i])) / (2 * M_PI * gsl_sf_bessel_I0(kappa[i]));
	}
	return result;
}

double DPVonMises::lnhyp(double x) {
	//return log(1 + x/2 + sqrt(x*x + 4)/2);
	return log(x*x + 1);
}

double DPVonMises::hyp(double x) {
	return sqrt(x*x + 1)-1;
}

void DPVonMises::print() {
	std::cout << "[";
	for (int i = 0; i < count; ++i) {
		std::cout << "(" << kappa[i] << ", " << mu[i] << ", " << weight[i] << "), ";
	}
	std::cout << std::endl;
}

int DPVonMises::moment_f(const gsl_vector* x, void* params, gsl_vector* f) {
	RightSide* rs = (RightSide*) params;
	int c = rs->count*2/3;

	for (int i = 0; i < rs->count; ++i) {
		double y_re = 0;
		double y_im = 0;

		for (int j = 0; j < c; ++j) {
			double x_kappa  = lnhyp(gsl_vector_get(x, 3*j));
			double x_mu     = gsl_vector_get(x, 3*j + 1);
			double x_weight = hyp(gsl_vector_get(x, 3*j + 2));

			double foo = gsl_sf_bessel_In(i+1, x_kappa) / gsl_sf_bessel_I0(x_kappa);

			y_re += x_weight * foo * cos((i+1)*x_mu);
			y_im += x_weight * foo * sin((i+1)*x_mu);
		}

		gsl_vector_set (f, 2*i    , y_re - rs->moment_re[i]);
		gsl_vector_set (f, 2*i + 1, y_im - rs->moment_im[i]);

	}

	return GSL_SUCCESS;
}
