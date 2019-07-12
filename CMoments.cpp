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
	moment1_re(),
	moment1_im()
{
}

CMoments::RightSide::RightSide(double moment1_re_, double moment1_im_) :
	moment1_re(moment1_re_),
	moment1_im(moment1_im_)
{
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

CMoments::DensityParams::DensityParams() :
	kappa(),
	mu()
{
}

CMoments::DensityParams::DensityParams(double kappa_, double mu_) :
	kappa(kappa_),
	mu(mu_)
{
}

CMoments::DensityParams::DensityParams(RightSide& rs) :
	kappa(),
	mu()
{
	int status;
	int iter = 0;

	const size_t n = 2;
	gsl_multiroot_function f = {&CMoments::moment_f, n, &rs};

	gsl_vector *x = gsl_vector_alloc(n);
	gsl_vector_set(x, 0, 1.86435);
	gsl_vector_set(x, 1, 1.35135);

	gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_hybrids, 2);
	gsl_multiroot_fsolver_set (s, &f, x);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);

		std::cout << "iter "<< iter <<": kappa = "<< gsl_vector_get (s->x, 0) <<", mu = "<< gsl_vector_get (s->x, 1) << std::endl;

		if (status) {   /* check if solver is stuck */
			break;
		}

		status = gsl_multiroot_test_residual (s->f, 1e-7);
	}	while (status == GSL_CONTINUE && iter < 1000);

	printf ("status = %s\n", gsl_strerror (status));
	kappa = gsl_vector_get (s->x, 0);
	mu = gsl_vector_get (s->x, 0) + M_PI_2;

	gsl_multiroot_fsolver_free(s);
}

double CMoments::DensityParams::density_at(double phase) {
	return exp(kappa * cos(phase - mu)) / (2 * M_PI * gsl_sf_bessel_I0(kappa));
}

CMoments::MomentEstimator::MomentEstimator() :
	sum1_re(),
	sum1_im(),
	count()
{
}

void CMoments::MomentEstimator::add_point(double phase) {
	sum1_re += cos(phase);
	sum1_im += sin(phase);
	count++;
}

CMoments::RightSide CMoments::MomentEstimator::estimate_moments() {
	return RightSide(sum1_re / count, sum1_im / count);
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

	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);

	double foo = gsl_sf_bessel_I1(x0) / gsl_sf_bessel_I0(x0);

	const double y0 = foo * cos(x1) - rs->moment1_re;
	const double y1 = foo * sin(x1) - rs->moment1_im;

	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);

	return GSL_SUCCESS;
}

void CMoments::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	RightSide pos_rs = pos_estimator.estimate_moments();
	RightSide neg_rs = neg_estimator.estimate_moments();

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
		std::cout << "positive: k="<< pos_density.kappa <<" m="<< pos_density.mu << std::endl;
		std::cout << "negative: k="<< neg_density.kappa <<" m="<< neg_density.mu << std::endl;
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
	fwrite(&pos_density.kappa, sizeof(double), 1, file);
	fwrite(&pos_density.mu, sizeof(double), 1, file);
	fwrite(&neg_density.kappa, sizeof(double), 1, file);
	fwrite(&neg_density.mu, sizeof(double), 1, file);
	return 0;
}

int CMoments::load(FILE* file)
{
	fread(&pos_density.kappa, sizeof(double), 1, file);
	fread(&pos_density.mu, sizeof(double), 1, file);
	fread(&neg_density.kappa, sizeof(double), 1, file);
	fread(&neg_density.mu, sizeof(double), 1, file);
	return 0;
}


int CMoments::exportToArray(double* array, int maxLen)
{
	int pos = 0;
	array[pos++] = type;
	array[pos++] = pos_density.kappa;
	array[pos++] = pos_density.mu;
	array[pos++] = neg_density.kappa;
	array[pos++] = neg_density.mu;
	return pos;
}

int CMoments::importFromArray(double* array, int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_NONE) std::cerr << "Error loading the model, type mismatch." << std::endl;
	pos_density.kappa = array[pos++];
	pos_density.mu = array[pos++];
	neg_density.kappa = array[pos++];
	neg_density.mu = array[pos++];
	return pos;
}
