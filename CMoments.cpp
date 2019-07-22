#include <iostream>
#include <fstream>
#include "CMoments.h"

#include "moments/right_side.h"

using namespace std;

CMoments::CMoments(int idd) :
	pos_estimator(),
	neg_estimator(),
	pos_density(),
	neg_density()
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

void CMoments::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	RightSide pos_rs(pos_estimator);
	RightSide neg_rs(neg_estimator);

	pos_density.calculate(pos_rs);
	neg_density.calculate(neg_rs);

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
	//float nd = neg_density.density_at(phase);

	//return pd / (pd + nd);
	return pd;
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
	fwrite(&pos_density.kappa[0], sizeof(double), pos_density.count, file);
	fwrite(&pos_density.mu[0], sizeof(double), pos_density.count, file);
	fwrite(&pos_density.weight[0], sizeof(double), pos_density.count, file);
	fwrite(&neg_density.count, sizeof(int), 1, file);
	fwrite(&neg_density.kappa[0], sizeof(double), neg_density.count, file);
	fwrite(&neg_density.mu[0], sizeof(double), neg_density.count, file);
	fwrite(&neg_density.weight[0], sizeof(double), neg_density.count, file);
	return 0;
}

int CMoments::load(FILE* file)
{
	int cnt;
	fread(&cnt, sizeof(int), 1, file);
	pos_density.reset(cnt);
	fread(&pos_density.kappa[0], sizeof(double), cnt, file);
	fread(&pos_density.mu[0], sizeof(double), cnt, file);
	fread(&pos_density.weight[0], sizeof(double), cnt, file);
	fread(&cnt, sizeof(int), 1, file);
	neg_density.reset(cnt);
	fread(&neg_density.kappa[0], sizeof(double), cnt, file);
	fread(&neg_density.mu[0], sizeof(double), cnt, file);
	fread(&neg_density.weight[0], sizeof(double), cnt, file);
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
	pos_density.reset(cnt);
	for (int i = 0; i <= pos_density.count; ++i) {
		pos_density.kappa[i] = array[pos++];
		pos_density.mu[i] = array[pos++];
		pos_density.weight[i] = array[pos++];
	}
	cnt = array[pos++];
	neg_density.reset(cnt);
	for (int i = 0; i <= pos_density.count; ++i) {
		neg_density.kappa[i] = array[pos++];
		neg_density.mu[i] = array[pos++];
		neg_density.weight[i] = array[pos++];
	}
	return pos;
}
