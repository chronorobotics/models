#include <iostream>
#include <fstream>
#include <ctime>
#include "CExpectation.h"

CExpectation::CExpectation(int idd) :
	models(),
	means(),
	class_model(0, 0)
{
	throw "We don't do this here";
}

CExpectation::CExpectation(int idd, int dimension_) :
	dimension(dimension_),
	models(),
	means(),
	class_model(5, dimension_)
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

CExpectation::TimeSample::TimeSample(long t_, std::vector<bool> v_) :
	t(t_),
	v(v_)
{
}

// adds new state observations at given times
int CExpectation::add(uint32_t time, float state)
{
	throw "We don't do this here";
}

int CExpectation::add_v(uint32_t time, std::vector<bool> state)
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
	class_model.print();
	means = class_model.get_means();
	std::cout << numSamples << " " << means.size() << " " << models.size() << std::endl;
	for (int i = 0; i < numSamples; ++i) {
		std::vector<double> alpha = class_model.get_alpha_at(sampleArray[i].v);
		for (int j = 0; j < alpha.size(); ++j) {
			models[j].add_time(sampleArray[i].t, alpha[j]);
		}
	}

	for (int i = 0; i < models.size(); ++i) {
		models[i].train();
	}

	//class_model.print();
	std::cout << "Trained" << std::endl;
}

/*text representation of the fremen model*/
void CExpectation::print(bool verbose)
{
	/*for (int i = 0; i < means.size(); ++i) {
		std::cout << "Mean (";
		for (int j = 0; j < dimension; ++j) {
			if (j) {
				std::cout << ", ";
			}
			std::cout << means[i][j];
		}
		std::cout << "): ";
		models[i].print();
		std::cout << std::endl;
	}*/
	class_model.print();
}

std::vector<float> CExpectation::estimate_v(uint32_t time)
{
	std::vector<float> s;
	s.resize(dimension, 0);
	double t = 0;
	for (int i = 0; i < means.size(); ++i) {
		double density = models[i].get_density_at(time);
		for (int j = 0; j < dimension; ++j) {
			s[j] += means[i][j] * density;
		}
		t += density;
	}

	for (int i = 0; i < dimension; ++i) {
		s[i] /= t;
	}

	return s;
}

float CExpectation::estimate(uint32_t time)
{
	throw "We don't do this here.";
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
	fwrite(&dimension, sizeof(int), 1, file);
	for (int i = 0; i < models.size(); ++i) {
		for (int j = 0; j < dimension; ++j) {
			fwrite(&(means[i][j]), sizeof(double), 1, file);
		}
		models[i].save(file, lossy);
	}
	return 0;
}

int CExpectation::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	fread(&dimension, sizeof(int), 1, file);
	means.resize(size);
	models.resize(size);
	for (int i = 0; i < models.size(); ++i) {
		means[i].resize(dimension);
		for (int j = 0; j < dimension; ++j) {
			fread(&(means[i][j]), sizeof(double), 1, file);
		}
		models[i].load(file);
	}
	return 0;
}

int CExpectation::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = TT_EXPECTATION;
	array[pos++] = models.size();
	array[pos++] = dimension;
	for (int i = 0; i < models.size(); ++i) {
		for (int j = 0; j < dimension; ++j) {
			array[pos++] = means[i][j];
		}
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
	dimension = array[pos++];
	means.resize(size);
	models.resize(size);
	for (int i = 0; i < models.size(); ++i) {
		means[i].resize(dimension);
		for (int j = 0; j < dimension; ++j) {
			means[i][j] = array[pos++];
		}
		models[i].importFromArray(array, len, pos);
	}
	return pos;
}

