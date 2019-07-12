#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>
#include "CMoments.h"

using namespace std;

CMoments::CMoments(int idd)
{
	id=idd;
	measurements = 0;
	type = TT_MOMENTS;

	numSamples = 0;
}

void CMoments::init(int iMaxPeriod,int elements,int numClasses)
{
	maxPeriod = iMaxPeriod;
	numElements = 1;
}

CMoments::~CMoments()
{
}

CMoments::RightSide::RightSide() :
	moment_re(),
	moment_im(),
	count(moment_count)
{
	moment_re = new double[moment_count];
	moment_im = new double[moment_count];
}

CMoments::RightSide::RightSide(const CMoments::MomentEstimator& me) :
	moment_re(),
	moment_im(),
	count(moment_count)
{
	moment_re = new double[moment_count];
	moment_im = new double[moment_count];

	for (int i = 0; i < moment_count; ++i) {
		moment_re[i] = me.sum_re[i] / me.count;
		moment_im[i] = me.sum_im[i] / me.count;
	}
}

CMoments::RightSide::~RightSide() {
	delete[] moment_re;
	delete[] moment_im;
}

CMoments::TimeSample::TimeSample() :
	t(),
	v()
{
}

CMoments::TimeSample::TimeSample(long t_, float v_) :
	t(t_),
	v(v_)
{
}

CMoments::DensityParams::DensityParams(int count_) :
	kappa(),
	mu(),
	weight(),
	count(count_)
{
}

CMoments::DensityParams::~DensityParams() {
	if (count) {
		delete[] kappa;
		delete[] mu;
		delete[] weight;
	}
}

CMoments::DensityParams::DensityParams(RightSide& rs) :
	kappa(),
	mu(),
	weight(),
	count()
{
	count = moment_count*2/3;
	kappa = new double[count];
	mu = new double[count];
	weight = new double[count];

	int status;
	int iter = 0;

	const size_t n = moment_count*2;
	std::cout << n << " " << count << std::endl;
	gsl_multiroot_function f = {&CMoments::moment_f, n, &rs};

	gsl_vector *x = gsl_vector_alloc(n);
	for (int i = n-1; i >= 0; --i) {
		gsl_vector_set(x, i, 1);
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

		double abs = 0;
		for (int i = n-1; i >= 0; --i) {
			double foo = gsl_vector_get(s->f, i);
			abs += foo*foo;
		}

		std::cout << "iter "<< iter <<": abs = "<< sqrt(abs) << std::endl;

		if (status) {
			break;
		}

		status = gsl_multiroot_test_residual (s->f, 0.07);
	}	while (status == GSL_CONTINUE && iter < 1000);

	printf ("status = %s\n", gsl_strerror (status));
	for (int i = 0; i < count; ++i) {
		kappa[i]  = gsl_vector_get(s->x, 3*i);
		mu[i]     = gsl_vector_get(s->x, 3*i + 1);
		weight[i] = gsl_vector_get(s->x, 3*i + 2);
	}

	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
}

double CMoments::DensityParams::density_at(double phase) {
	double result = 0;
	for (int i = 0; i < count; ++i) {
		result += weight[i] * exp(kappa[i] * cos(phase - mu[i])) / (2 * M_PI * gsl_sf_bessel_I0(kappa[i]));
	}
	return result;
}

void CMoments::DensityParams::print() {
	std::cout << "[";
	for (int i = 0; i < count; ++i) {
		std::cout << "(" << kappa[i] << ", " << mu[i] << ", " << weight[i] << "), ";
	}
	std::cout << std::endl;
}

CMoments::MomentEstimator::MomentEstimator() :
	sum_re(),
	sum_im(),
	count(0)
{
	sum_re = new double[moment_count];
	sum_im = new double[moment_count];
}

CMoments::MomentEstimator::~MomentEstimator() {
	delete[] sum_re;
	delete[] sum_im;
}

void CMoments::MomentEstimator::add_point(double phase) {
	for (int i = 1; i <= moment_count; ++i) {
		sum_re[i] += cos(i*phase);
		sum_im[i] += sin(i*phase);
	}
	count++;
}

double CMoments::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 86400.0f) / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

// adds new state observations at given times
int CMoments::add(uint32_t time,float state)
{
	sampleArray[numSamples] = TimeSample(time, state);
	numSamples++;

	float phase = time_to_phase(time);

	if (state > 0.5) {
		pos_estimator.add_point(phase);
	} else {
		neg_estimator.add_point(phase);
	}
	measurements++;
	return 0;
}

int CMoments::moment_f(const gsl_vector* x, void* params, gsl_vector* f) {
	RightSide* rs = (RightSide*) params;
	int c = rs->count*2/3;

	/*
		const double x0 = gsl_vector_get (x, 0);
		const double x1 = gsl_vector_get (x, 1);

		double foo = gsl_sf_bessel_I1(x0) / gsl_sf_bessel_I0(x0);

		const double y0 = foo * cos(x1) - rs->moment1_re;
		const double y1 = foo * sin(x1) - rs->moment1_im;

		gsl_vector_set (f, 0, y0);
		gsl_vector_set (f, 1, y1);
	*/

	for (int i = 0; i < rs->count; ++i) {
		double y_re = 0;
		double y_im = 0;

		for (int j = 0; j < c; ++j) {
			const double x_kappa  = gsl_vector_get(x, 3*j);
			const double x_mu     = gsl_vector_get(x, 3*j + 1);
			const double x_weight = gsl_vector_get(x, 3*j + 2);

			double foo = gsl_sf_bessel_In(i+1, x_kappa) / gsl_sf_bessel_I0(x_kappa);

			y_re += x_weight * foo * cos((i+1)*x_mu);
			y_im += x_weight * foo * sin((i+1)*x_mu);
		}

		gsl_vector_set (f, 2*i    , y_re - rs->moment_re[i]);
		gsl_vector_set (f, 2*i + 1, y_im - rs->moment_im[i]);
	}

	return GSL_SUCCESS;
}

void CMoments::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	RightSide pos_rs(pos_estimator);
	RightSide neg_rs(neg_estimator);

	pos_density = DensityParams(pos_rs);
	neg_density = DensityParams(neg_rs);

	ofstream myfile0("0.txt");
	ofstream myfile1("1.txt");
	for (int i = 0; i < numSamples; ++i) {
		myfile0 << sampleArray[i].v << std::endl;
		myfile1 << predict(sampleArray[i].t) << std::endl;
	}
	myfile0.close();
	myfile1.close();
}

/*text representation of the fremen model*/
void CMoments::print(bool verbose)
{
	std::cout << "Model " << id << " Size: " << measurements << " " << std::endl;
	if (verbose) {
		std::cout << "positive: ";
		pos_density.print();
		std::cout << std::endl << "negative: ";
		neg_density.print();
		std::cout << std::endl;
	}
}

float CMoments::estimate(uint32_t time)
{
	float phase = time_to_phase(time);

	float pd = pos_density.density_at(phase);
	float nd = neg_density.density_at(phase);

	return pd / (pd + nd);
	//return pd;
}

float CMoments::predict(uint32_t time)
{
	return estimate(time);
}

int CMoments::save(const char* name, bool lossy)
{
	FILE* file = fopen(name,"w");
	save(file);
	fclose(file);
	return 0;
}

int CMoments::load(const char* name)
{
	FILE* file = fopen(name,"r");
	load(file);
	fclose(file);
	return 0;
}


int CMoments::save(FILE* file,bool lossy)
{
	fwrite(&pos_density.count, sizeof(int), 1, file);
	fwrite(pos_density.kappa, sizeof(double), pos_density.count, file);
	fwrite(pos_density.mu, sizeof(double), pos_density.count, file);
	fwrite(pos_density.weight, sizeof(double), pos_density.count, file);
	fwrite(&neg_density.count, sizeof(int), 1, file);
	fwrite(neg_density.kappa, sizeof(double), neg_density.count, file);
	fwrite(neg_density.mu, sizeof(double), neg_density.count, file);
	fwrite(neg_density.weight, sizeof(double), neg_density.count, file);
	return 0;
}

int CMoments::load(FILE* file)
{
	int cnt;
	fread(&cnt, sizeof(int), 1, file);
	pos_density = DensityParams(cnt);
	fread(pos_density.kappa, sizeof(double), cnt, file);
	fread(pos_density.mu, sizeof(double), cnt, file);
	fread(pos_density.weight, sizeof(double), cnt, file);
	fread(&cnt, sizeof(int), 1, file);
	neg_density = DensityParams(cnt);
	fread(neg_density.kappa, sizeof(double), cnt, file);
	fread(neg_density.mu, sizeof(double), cnt, file);
	fread(neg_density.weight, sizeof(double), cnt, file);
	return 0;
}


int CMoments::exportToArray(double* array, int maxLen)
{
	int pos = 0;
	array[pos++] = type;
	array[pos++] = pos_density.count;
	for (int i = 0; i <= pos_density.count; ++i) {
		array[pos++] = pos_density.kappa[i];
		array[pos++] = pos_density.mu[i];
		array[pos++] = pos_density.weight[i];
	}
	array[pos++] = neg_density.count;
	for (int i = 0; i <= neg_density.count; ++i) {
		array[pos++] = neg_density.kappa[i];
		array[pos++] = neg_density.mu[i];
		array[pos++] = neg_density.weight[i];
	}
	return pos;
}

int CMoments::importFromArray(double* array, int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_NONE) std::cerr << "Error loading the model, type mismatch." << std::endl;
	int cnt;
	cnt = array[pos++];
	pos_density = DensityParams(cnt);
	for (int i = 0; i <= pos_density.count; ++i) {
		pos_density.kappa[i] = array[pos++];
		pos_density.mu[i] = array[pos++];
		pos_density.weight[i] = array[pos++];
	}
	cnt = array[pos++];
	neg_density = DensityParams(cnt);
	for (int i = 0; i <= pos_density.count; ++i) {
		neg_density.kappa[i] = array[pos++];
		neg_density.mu[i] = array[pos++];
		neg_density.weight[i] = array[pos++];
	}
	return pos;
}
