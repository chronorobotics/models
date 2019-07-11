#include <gsl/gsl_sf_bessel.h>
#include "CMoments.h"

using namespace std;

CMoments::CMoments(int idd)
{
	id=idd;
	measurements = 0;
	type = TT_MOMENTS;
	pos_sum1_re = 0;
	pos_sum1_im = 0;
	neg_sum1_re = 0;
	neg_sum1_im = 0;
	positives = 0;
	negatives = 0;
}

void CMoments::init(int iMaxPeriod,int elements,int numClasses)
{
	maxPeriod = iMaxPeriod;
	numElements = 1;
}

CMoments::~CMoments()
{
}

// adds new state observations at given times
int CMoments::add(uint32_t time,float state)
{
	float phase = time / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	phase *= M_PI;

	if (state > 0.5) {
		pos_sum1_re += cos(phase);
		pos_sum1_im += sin(phase);
		positives++;
	} else {
		neg_sum1_re += cos(phase);
		neg_sum1_im += sin(phase);
		negatives++;
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
	RightSide pos_rs, neg_rs;
	pos_rs.moment1_re = pos_sum1_re / positives;
	pos_rs.moment1_im = pos_sum1_im / positives;
	neg_rs.moment1_re = neg_sum1_re / negatives;
	neg_rs.moment1_im = neg_sum1_im / negatives;

	int status;
	int iter = 0;

	const size_t n = 2;
	gsl_multiroot_function f_pos = {&moment_f, n, &pos_rs};
	gsl_multiroot_function f_neg = {&moment_f, n, &neg_rs};

	gsl_vector *x = gsl_vector_alloc(n);
	gsl_vector_set(x, 0, 1.86435);
	gsl_vector_set(x, 1, 1.35135);

	gsl_multiroot_fsolver* s1 = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_hybrids, 2);
	gsl_multiroot_fsolver_set (s1, &f_pos, x);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s1);

		std::cout << "iter "<< iter <<": kappa = "<< gsl_vector_get (s1->x, 0) <<", mu = "<< gsl_vector_get (s1->x, 1) << std::endl;

		if (status) {   /* check if solver is stuck */
			break;
		}

		status = gsl_multiroot_test_residual (s1->f, 1e-7);
	}	while (status == GSL_CONTINUE && iter < 1000);

	printf ("status = %s\n", gsl_strerror (status));
	pos_kappa = gsl_vector_get (s1->x, 0);
	pos_mu = gsl_vector_get (s1->x, 0);

	gsl_multiroot_fsolver_free(s1);

	gsl_multiroot_fsolver* s2 = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_hybrids, 2);
	gsl_multiroot_fsolver_set (s2, &f_neg, x);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s2);

		std::cout << "iter "<< iter <<": kappa = "<< gsl_vector_get (s2->x, 0) <<", mu = "<< gsl_vector_get (s2->x, 1) << std::endl;

		if (status) {   /* check if solver is stuck */
			break;
		}

		status =gsl_multiroot_test_residual (s2->f, 1e-7);
	}	while (status == GSL_CONTINUE && iter < 1000);

	printf ("status = %s\n", gsl_strerror (status));
	neg_kappa = gsl_vector_get (s2->x, 0);
	neg_mu = gsl_vector_get (s2->x, 0);

	gsl_multiroot_fsolver_free(s2);
	gsl_vector_free(x);
}

/*text representation of the fremen model*/
void CMoments::print(bool verbose)
{
	std::cout << "Model " << id << " Size: " << measurements << " " << std::endl;
	if (verbose) {
		std::cout << "positive: k="<< pos_kappa <<" m="<< pos_mu << std::endl;
		std::cout << "negative: k="<< neg_kappa <<" m="<< neg_mu << std::endl;
	}
}

float CMoments::estimate(uint32_t time)
{
	float phase = time / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	phase *= M_PI;

	float pos_density = exp(pos_kappa * cos(phase - pos_mu)) / (2 * M_PI * gsl_sf_bessel_I0(pos_kappa));
	float neg_density = exp(neg_kappa * cos(phase - neg_mu)) / (2 * M_PI * gsl_sf_bessel_I0(neg_kappa));

	return pos_density / (pos_density + neg_density);
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
	fwrite(&pos_kappa, sizeof(double), 1, file);
	fwrite(&pos_mu, sizeof(double), 1, file);
	fwrite(&neg_kappa, sizeof(double), 1, file);
	fwrite(&neg_mu, sizeof(double), 1, file);
	return 0;
}

int CMoments::load(FILE* file)
{
	fread(&pos_kappa, sizeof(double), 1, file);
	fread(&pos_mu, sizeof(double), 1, file);
	fread(&neg_kappa, sizeof(double), 1, file);
	fread(&neg_mu, sizeof(double), 1, file);
	return 0;
}


int CMoments::exportToArray(double* array, int maxLen)
{
	int pos = 0;
	array[pos++] = type;
	array[pos++] = estimation;
	array[pos++] = id;
	array[pos++] = measurements;
	return pos;
}

int CMoments::importFromArray(double* array, int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_NONE) std::cerr << "Error loading the model, type mismatch." << std::endl;
	estimation = array[pos++];
	id = array[pos++];  
	measurements = array[pos++]; 
	return pos;
}
