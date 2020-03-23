#include <iostream>
#include <fstream>
#include <ctime>
#include "CAdaboost.h"
#include "CFrelement.h"

CAdaboost::CAdaboost(int idd) :
	adaboosts(),
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
		adaboosts.push_back(Adaboost());
	}

	srandom(time(0));
}

void CAdaboost::init(int iMaxPeriod, int elements, int numClasses)
{
	maxPeriod = iMaxPeriod;
	numElements = 1;
}

CAdaboost::~CAdaboost()
{
}


CAdaboost::TimeSample::TimeSample() :
	t(),
	v()
{
}

CAdaboost::TimeSample::TimeSample(long t_, float v_) :
	t(t_),
	v(v_)
{
}

// adds new state observations at given times
int CAdaboost::add(uint32_t time, float state)
{
	/*sampleArray[numSamples] = TimeSample(time, state);
	numSamples++;
	my_adaboost.add_time(time, state > 0.5);*/
	sampleArray[numSamples] = TimeSample(time, state);
	numSamples++;

	class_model.add_value(state);
	return 0;
}

/*not required in incremental version*/
void CAdaboost::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	class_model.train();
	means = class_model.get_means();
	std::cout << numSamples << " " << means.size() << " " << adaboosts.size() << std::endl;
	for (int i = 0; i < numSamples; ++i) {
		std::vector<double> alpha = class_model.get_alpha_at(sampleArray[i].v);
		double max = 0;
		int argmax = 0;
		for (int j = 0; j < alpha.size(); ++j) {
			if (alpha[j] > max) {
				max = alpha[j];
				argmax = j;
			}
		}
		for (int j = 0; j < alpha.size(); ++j) {
			adaboosts[j].add_time(sampleArray[i].t, j==argmax);
		}
	}

	CFrelement fremen(0);
	fremen.init(maxPeriod, id, 1);
	for (int i = 0; i < numSamples; ++i) {
		fremen.add(sampleArray[i].t, sampleArray[i].v);
	}
	fremen.update(id);
	for (int i = 0; i < id; ++i) {
		for (int j = 0; j < adaboosts.size(); ++j) {
			adaboosts[j].add_period(fremen.getPredictFrelements()[i].period);
		}
	}
	for (int j = 0; j < adaboosts.size(); ++j) {
		adaboosts[j].train();
	}

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
void CAdaboost::print(bool verbose)
{
	for (int j = 0; j < adaboosts.size(); ++j) {
		adaboosts[j].print();
	}
	//my_adaboost.print();
}

float CAdaboost::estimate(uint32_t time)
{
	//return (my_adaboost.classify(time)+1)/2;
	double s = 0;
	double t = 0;
	for (int i = 0; i < means.size(); ++i) {
		double density = (adaboosts[i].classify(time)+1)/2;
		s += means[i] * density;
		t += density;
	}

	return s/t;
}

float CAdaboost::predict(uint32_t time)
{
	return estimate(time);
}

int CAdaboost::save(const char* name, bool lossy)
{
	FILE* file = fopen(name,"w");
	save(file);
	fclose(file);
	return 0;
}

int CAdaboost::load(const char* name)
{
	FILE* file = fopen(name,"r");
	load(file);
	fclose(file);
	return 0;
}


int CAdaboost::save(FILE* file, bool lossy)
{
	//my_adaboost.save(file, lossy);
	int size = adaboosts.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < adaboosts.size(); ++i) {
		fwrite(&(means[i]), sizeof(double), 1, file);
		adaboosts[i].save(file, lossy);
	}
	return 0;
}

int CAdaboost::load(FILE* file)
{
	//my_adaboost.load(file);
	int size;
	fread(&size, sizeof(int), 1, file);
	means.resize(size);
	adaboosts.resize(size);
	for (int i = 0; i < adaboosts.size(); ++i) {
		fread(&(means[i]), sizeof(double), 1, file);
		adaboosts[i].load(file);
	}
	return 0;
}

int CAdaboost::exportToArray(double* array,int maxLen)
{
	/*int pos = 0;
	array[pos++] = type;
	my_adaboost.exportToArray(array, maxLen, pos);*/
	int pos = 0;
	array[pos++] = adaboosts.size();
	for (int i = 0; i < adaboosts.size(); ++i) {
		array[pos++] = means[i];
		adaboosts[i].exportToArray(array, maxLen, pos);
	}
	return pos;
}

int CAdaboost::importFromArray(double* array,int len)
{
	/*int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_ADABOOST) std::cerr << "Error loading the model, type mismatch." << std::endl;
	my_adaboost.importFromArray(array, len, pos);*/
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_ADABOOST) std::cerr << "Error loading the model, type mismatch." << std::endl;
	int size = array[pos++];
	means.resize(size);
	adaboosts.resize(size);
	for (int i = 0; i < adaboosts.size(); ++i) {
		means[i] = array[pos++];
		adaboosts[i].importFromArray(array, len, pos);
	}
	return pos;
}

