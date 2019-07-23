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
	int status;
	int tries = 0;
	srandom(time(0));

	do {
		tries++;
		int iter = 0;
		const size_t n = ep.right_side.size();
		gsl_multiroot_function f = {&DPVonMises::moment_f, n, &ep};

		gsl_vector* x = gsl_vector_alloc(n);
		for (int i = 0; i < n; ++i) {
			gsl_vector_set(x, i, float(rand()) / RAND_MAX * 2);
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

			/*std::cout << "iter " << iter << std::endl;
			std::cout << "    residuum:";
			for (int i = 0; i < n; ++i) {
				std::cout << " " << gsl_vector_get(s->f, i);
			}
			std::cout << std::endl << "    value:";
			for (int i = 0; i < n; ++i) {
				std::cout << " " << gsl_vector_get(s->x, i);
			}
			std::cout << std::endl;*/

			if (status) {
				break;
			}

			status = gsl_multiroot_test_residual (s->f, 1E-7);
		}	while (status == GSL_CONTINUE && iter < 1000);

		std::cout << "status = " << gsl_strerror (status) << ", tries = " << tries << std::endl;
		for (int i = 0; i < count; ++i) {
			kappa[i]  = lnhyp(gsl_vector_get(s->x, 3*i));
			mu[i]     = gsl_vector_get(s->x, 3*i + 1);
			weight[i] = hyp(gsl_vector_get(s->x, 3*i + 2));
		}

		gsl_multiroot_fsolver_free(s);
		gsl_vector_free(x);
	} while (status != GSL_SUCCESS && tries < 100);

	if (status != GSL_SUCCESS) {
		std::cout << "ERROR: Solution not found!" << std::endl;
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
	EquationParams* ep = (EquationParams*) params;

	for (int i = ep->moment_count - 1; i >= 0; --i) {
		double y_re = 0;
		double y_im = 0;

		for (int j = ep->cluster_count - 1; j >= 0; --j) {
			double x_kappa  = lnhyp(gsl_vector_get(x, 3*j));
			double x_mu     = gsl_vector_get(x, 3*j + 1);
			double x_weight = hyp(gsl_vector_get(x, 3*j + 2));

			if (isnan(x_kappa) || isnan(x_mu) || isnan(x_weight)) {
				return GSL_EDOM;
			}

			double foo = gsl_sf_bessel_In(i+1, x_kappa) / gsl_sf_bessel_I0(x_kappa);

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
	moment_count(moment_count_)
{

}
