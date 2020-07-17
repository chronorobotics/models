#include <math.h>
#include <ctime>

#include "../CMoments.h"
#include "me_circular.h"
#include "dp_sqdist.h"

DPSqdist::DPSqdist(CMoments* parent_, int count_) :
	DensityParams(parent_, count_),
	xx(),
	yy(),
	weight(),
	estimator()
{
	xx.resize(count);
	yy.resize(count);
	weight.resize(count);
}

DPSqdist::~DPSqdist() {

}

MomentEstimator* DPSqdist::get_moment_estimator() {
	if (!estimator.get()) {
		estimator = std::shared_ptr<MECircular>(new MECircular(ceil(float(count)*3/2)));
	}
	return estimator.get();
}

void DPSqdist::reset(int count_) {
	count = count_;
	xx.resize(count);
	yy.resize(count);
	weight.resize(count);
}

void DPSqdist::calculate()
{
	EquationParams ep(estimator->get_moments(), count, estimator->get_moment_count());
	int status;
	int tries = 0;
	srandom(time(0));

	if (count == 1) {
		xx[0] = ep.right_side[0];
		yy[0] = ep.right_side[1];
		weight[0] = 1;
		return;
	}

	do {
		tries++;
		int iter = 0;
		const size_t n = count*3;
		gsl_multiroot_function f = {&DPSqdist::moment_f, n, &ep};

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

			status = gsl_multiroot_test_residual (s->f, 0.1);
		}	while (status == GSL_CONTINUE && iter < 1000);

		std::cout << "\33[2K\rstatus = " << gsl_strerror (status) << ", tries = " << tries << std::endl;
		double sum_w = 0;
		bool negative_weight = false;
		for (int i = 0; i < count; ++i) {
			xx[i]     = gsl_vector_get(s->x, 3*i);
			yy[i]     = gsl_vector_get(s->x, 3*i + 1);
			weight[i] = gsl_vector_get(s->x, 3*i + 2);
			sum_w += weight[i];
			if (weight[i] < 0) {
				negative_weight = true;
			}
		}

		/*for (int i = 0; i < count; ++i) {
			weight[i] /= sum_w;
		}*/

		if (status == GSL_SUCCESS) {
			for (int i = 0; i < n; ++i) {
				std::cout << (gsl_vector_get(s->f, 0) + ep.right_side[i]) << ", ";
				if (i % 2) {
					std::cout << std::endl;
				}
			}
			std::cout << std::endl;

			for (int i = 0; i < n; ++i) {
				std::cout << ep.right_side[i] << ", ";
				if (i % 2) {
					std::cout << std::endl;
				}
			}
			std::cout << std::endl;
		}

		gsl_multiroot_fsolver_free(s);
		gsl_vector_free(x);
		if (sum_w > 1.1 || sum_w < 0.9 || negative_weight) {
			status = GSL_EDOM;
		}
	} while (status != GSL_SUCCESS);

	std::cout << std::endl;

	if (status != GSL_SUCCESS) {
		std::cout << "ERROR: Solution not found!" << std::endl;
	}
}

double DPSqdist::density_at(uint32_t time) {
	double phase = MECircular::time_to_phase(time);
	double result = 0;
	for (int i = 0; i < count; ++i) {
		double x = cos(phase) - xx[i];
		double y = sin(phase) - yy[i];
		result += weight[i] * (1 - (xx[i]*xx[i] + yy[i]*yy[i])) / (2*M_PI * (x*x + y*y));
	}
	return result;
}

void DPSqdist::print() {
	std::cout << "[";
	for (int i = 0; i < count; ++i) {
		std::cout << "(" << xx[i] << ", " << yy[i] << ", " << weight[i] << "), ";
	}
	std::cout << std::endl;
}

void DPSqdist::power(double& re, double& im, int n) {
	double r = 1;
	double i = 0;
	double r2, i2;

	while (n) {
		if (n % 2) {
			r2 = r*re - i*im;
			i2 = r*im + i*re;
			r = r2;
			i = i2;
		}

		r2 = re*re - im*im;
		i2 = 2*re*im;
		re = r2;
		im = i2;

		n >>= 1;
	}

	re = r;
	im = i;
}

int DPSqdist::moment_f(const gsl_vector* x, void* params, gsl_vector* f) {
	EquationParams* ep = (EquationParams*) params;

	for (int i = ep->moment_count - 1; i >= 0; --i) {
		double y_re = 0;
		double y_im = 0;
		double sum_w = 0;

		for (int j = ep->cluster_count - 1; j >= 0; --j) {
			double x_xx     = gsl_vector_get(x, 3*j);
			double x_yy     = gsl_vector_get(x, 3*j + 1);
			double x_weight = gsl_vector_get(x, 3*j + 2);
			sum_w += x_weight;

			if (isnan(x_xx) || isnan(x_yy) || isnan(x_weight)) {
				return GSL_EDOM;
			}

			power(x_xx, x_yy, i+1);

			y_re += x_weight * x_xx;
			y_im += x_weight * x_yy;
		}

		gsl_vector_set (f, 2*i    , y_re - ep->right_side[2*i]);
		if (2*i + 1 < ep->cluster_count*3) {
			gsl_vector_set (f, 2*i + 1, y_im - ep->right_side[2*i + 1]);
		}

	}

	return GSL_SUCCESS;
}

int DPSqdist::save(FILE* file, bool lossy)
{
	Distribution type = SQDIST;
	fwrite(&type, sizeof(Distribution), 1, file);

	fwrite(&count, sizeof(int), 1, file);
	fwrite(&xx[0], sizeof(double), count, file);
	fwrite(&yy[0], sizeof(double), count, file);
	fwrite(&weight[0], sizeof(double), count, file);
	return 0;
}

int DPSqdist::load(FILE* file)
{
	int cnt;
	fread(&cnt, sizeof(int), 1, file);
	reset(cnt);
	fread(&xx[0], sizeof(double), cnt, file);
	fread(&yy[0], sizeof(double), cnt, file);
	fread(&weight[0], sizeof(double), cnt, file);
	return 0;
}


void DPSqdist::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = (double) SQDIST;
	array[pos++] = count;
	for (int i = 0; i <= count; ++i) {
		array[pos++] = xx[i];
		array[pos++] = yy[i];
		array[pos++] = weight[i];
	}
}

void DPSqdist::importFromArray(double* array, int len, int& pos)
{
	int cnt;
	cnt = array[pos++];
	reset(cnt);
	for (int i = 0; i <= count; ++i) {
		xx[i] = array[pos++];
		yy[i] = array[pos++];
		weight[i] = array[pos++];
	}
}

DPSqdist::EquationParams::EquationParams(std::vector<double> right_side_, int cluster_count_, int moment_count_) :
	right_side(right_side_),
	cluster_count(cluster_count_),
	moment_count(moment_count_)
{

}
