#include <gsl/gsl_sf_bessel.h>
#include <ctime>
#include "../CMoments.h"

#include "dp_uniform.h"
#include "me_circular.h"

DPUniform::DPUniform(CMoments* parent_, int count_) :
	DensityParams(parent_, count_),
	start(),
	end(),
	weight(),
	estimator()
{
	start.resize(count);
	end.resize(count);
	weight.resize(count);
}

DPUniform::~DPUniform() {

}

MomentEstimator* DPUniform::get_moment_estimator() {
	if (!estimator.get()) {
		estimator = std::shared_ptr<MECircular>(new MECircular(ceil(float(count)*3/2)));
	}
	return estimator.get();
}

void DPUniform::reset(int count_) {
	count = count_;
	start.resize(count);
	end.resize(count);
	weight.resize(count);
}

void DPUniform::calculate()
{
	EquationParams ep(estimator->get_moments(), count, estimator->get_moment_count());
	int status;
	int tries = 0;
	srandom(time(0));

	do {
		tries++;
		int iter = 0;
		const size_t n = count*3;
		gsl_multiroot_function f = {&DPUniform::moment_f, n, &ep};

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

			status = gsl_multiroot_test_residual (s->f, 0.01);
		}	while (status == GSL_CONTINUE && iter < 1000);

		std::cout << "\33[2K\rstatus = " << gsl_strerror (status) << ", tries = " << tries << std::flush;
		for (int i = 0; i < count; ++i) {
			start[i] = normalise_angle(gsl_vector_get(s->x, 3*i));
			end[i]  = normalise_angle(gsl_vector_get(s->x, 3*i + 1));
			weight[i] = gsl_vector_get(s->x, 3*i + 2);
		}

		gsl_multiroot_fsolver_free(s);
		gsl_vector_free(x);
	} while (status != GSL_SUCCESS /*&& tries < 100*/);

	std::cout << std::endl;

	if (status != GSL_SUCCESS) {
		std::cout << "ERROR: Solution not found!" << std::endl;
	}
}

double DPUniform::density_at(uint32_t time) {
	double phase = MECircular::time_to_phase(time);
	double result = 0;
	for (int i = 0; i < count; ++i) {
		if ((start[i] < phase && phase < end[i]) || ((end[i] < start[i]) && (phase > start[i] || phase < end[i]))) {
			result += weight[i] / arc_length(start[i], end[i]);
		}
	}
	return result;
}

void DPUniform::print() {
	std::cout << "[";
	for (int i = 0; i < count; ++i) {
		std::cout << "(" << start[i] << ", " << end[i] << ", " << weight[i] << "), ";
	}
	std::cout << std::endl;
}

double DPUniform::normalise_angle(double x) {
	x = fmod(x, 2*M_PI);

	if (x > M_PI) {
		x -= 2*M_PI;
	} else if (x < -M_PI) {
		x += 2*M_PI;
	}

	return x;
}

double DPUniform::arc_length(double a, double b) {
	double l = b - a;
	if (l < 0) {
		l += 2*M_PI;
	}
	return l;
}

int DPUniform::moment_f(const gsl_vector* x, void* params, gsl_vector* f) {
	EquationParams* ep = (EquationParams*) params;

	for (int i = ep->moment_count - 1; i >= 0; --i) {
		double y_re = 0;
		double y_im = 0;

		for (int j = ep->cluster_count - 1; j >= 0; --j) {
			double x_start  = normalise_angle(gsl_vector_get(x, 3*j));
			double x_end    = normalise_angle(gsl_vector_get(x, 3*j + 1));
			double x_weight = gsl_vector_get(x, 3*j + 2);

			if (isnan(x_start) || isnan(x_end) || isnan(x_weight)) {
				return GSL_EDOM;
			}

			double foo = x_weight / (arc_length(x_start, x_end) * (i+1));

			y_re += x_weight * foo * (sin((i+1)*x_end) - sin((i+1)*x_start));
			y_im += x_weight * foo * (cos((i+1)*x_start) - cos((i+1)*x_end));
		}

		gsl_vector_set (f, 2*i    , y_re - ep->right_side[2*i]);
		if (2*i + 1 < ep->cluster_count*3) {
			gsl_vector_set (f, 2*i + 1, y_im - ep->right_side[2*i + 1]);
		}

	}

	return GSL_SUCCESS;
}

int DPUniform::save(FILE* file, bool lossy)
{
	Distribution type = UNIFORM;
	fwrite(&type, sizeof(Distribution), 1, file);

	fwrite(&count, sizeof(int), 1, file);
	fwrite(&start[0], sizeof(double), count, file);
	fwrite(&end[0], sizeof(double), count, file);
	fwrite(&weight[0], sizeof(double), count, file);
	return 0;
}

int DPUniform::load(FILE* file)
{
	int cnt;
	fread(&cnt, sizeof(int), 1, file);
	reset(cnt);
	fread(&start[0], sizeof(double), cnt, file);
	fread(&end[0], sizeof(double), cnt, file);
	fread(&weight[0], sizeof(double), cnt, file);
	return 0;
}


void DPUniform::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = (double) UNIFORM;
	array[pos++] = count;
	for (int i = 0; i <= count; ++i) {
		array[pos++] = start[i];
		array[pos++] = end[i];
		array[pos++] = weight[i];
	}
}

void DPUniform::importFromArray(double* array, int len, int& pos)
{
	int cnt;
	cnt = array[pos++];
	reset(cnt);
	for (int i = 0; i <= count; ++i) {
		start[i] = array[pos++];
		end[i] = array[pos++];
		weight[i] = array[pos++];
	}
}

DPUniform::EquationParams::EquationParams(std::vector<double> right_side_, int cluster_count_, int moment_count_) :
	right_side(right_side_),
	cluster_count(cluster_count_),
	moment_count(moment_count_)
{

}

