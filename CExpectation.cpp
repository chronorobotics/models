#include <iostream>
#include <fstream>
#include <ctime>
#include "CExpectation.h"

CExpectation::CExpectation(int idd) :
	models(),
	means(),
	class_model(7)
{
	id=idd;
	firstTime = -1;
	lastTime = -1;
	measurements = 0;
	maxPeriod = 0;
	numElements = 0;
	type = TT_EXPECTATION;
	numSamples = 0;

	for (int i = class_model.get_cluster_count(); i; --i) {
		models.push_back(EMSqdist(idd));
	}

	srandom(time(0));
}

void CExpectation::init(int iMaxPeriod, int elements, int numClasses)
{
	maxPeriod = iMaxPeriod;
	numElements = 1;
}

CExpectation::~CExpectation()
{
}


CExpectation::TimeSample::TimeSample() :
	t(),
	v()
{
}

CExpectation::TimeSample::TimeSample(long t_, float v_) :
	t(t_),
	v(v_)
{
}

// adds new state observations at given times
int CExpectation::add(uint32_t time, float state)
{
	sampleArray[numSamples] = TimeSample(time, state);
	numSamples++;

	class_model.add_value(state);
	return 0;
}

/*not required in incremental version*/
void CExpectation::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	class_model.train();
	means = class_model.get_means();
	std::cout << numSamples << " " << means.size() << " " << models.size() << std::endl;
	for (int i = 0; i < numSamples; ++i) {
		std::vector<double> alpha = class_model.get_alpha_at(sampleArray[i].v);
		int index = 0;
		double max = 0;
		for (int j = 0; j < alpha.size(); ++j) {
			if (alpha[j] > max) {
				max = alpha[j];
				index = j;
			}
		}
		models[index].add_time(sampleArray[i].t, 1);
	}

	for (int i = 0; i < models.size(); ++i) {
		models[i].train();
	}

	class_model.print();

	ofstream myfile0("0.txt");
	ofstream myfile1("1.txt");
	//ofstream myfile2("2.txt");
	for (int i = 0; i < numSamples; ++i) {
		myfile0 << sampleArray[i].v << std::endl;
		myfile1 << estimate(sampleArray[i].t) << std::endl;
		//myfile1 << positive.get_density_at(sampleArray[i].t)/negatives << std::endl;
		//myfile2 << negative.get_density_at(sampleArray[i].t)/positives*100 << std::endl;
	}
	myfile0.close();
	myfile1.close();
	//myfile2.close();
}

/*text representation of the fremen model*/
void CExpectation::print(bool verbose)
{
	for (int i = 0; i < means.size(); ++i) {
		std::cout << "Mean " << means[i] << ": ";
		models[i].print();
		std::cout << std::endl;
	}
}

float CExpectation::estimate(uint32_t time)
{
	double s = 0;
	double t = 0;
	for (int i = 0; i < means.size(); ++i) {
		double density = models[i].get_density_at(time);
		s += means[i] * density;
		t += density;
	}

	return s/t;
}

float CExpectation::predict(uint32_t time)
{
	return estimate(time);
}

int CExpectation::save(const char* name, bool lossy)
{
	FILE* file = fopen(name,"w");
	save(file);
	fclose(file);
	return 0;
}

int CExpectation::load(const char* name)
{
	FILE* file = fopen(name,"r");
	load(file);
	fclose(file);
	return 0;
}


int CExpectation::save(FILE* file, bool lossy)
{
	int size = models.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < models.size(); ++i) {
		fwrite(&(means[i]), sizeof(double), 1, file);
		models[i].save(file, lossy);
	}
	return 0;
}

int CExpectation::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	means.resize(size);
	models.resize(size);
	for (int i = 0; i < models.size(); ++i) {
		fread(&(means[i]), sizeof(double), 1, file);
		models[i].load(file);
	}
	return 0;
}

int CExpectation::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = models.size();
	for (int i = 0; i < models.size(); ++i) {
		array[pos++] = means[i];
		models[i].exportToArray(array, maxLen, pos);
	}
	return pos;
}

int CExpectation::importFromArray(double* array,int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_EXPECTATION) std::cerr << "Error loading the model, type mismatch." << std::endl;
	int size = array[pos++];
	means.resize(size);
	models.resize(size);
	for (int i = 0; i < models.size(); ++i) {
		means[i] = array[pos++];
		models[i].importFromArray(array, len, pos);
	}
	update(0);
	return pos;
}

