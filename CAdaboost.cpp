#include <iostream>
#include <fstream>
#include <ctime>
#include "CAdaboost.h"
#include "CFrelement.h"

CAdaboost::CAdaboost(int idd) :
	my_adaboost()
{
	id=idd;
	firstTime = -1;
	lastTime = -1;
	measurements = 0;
	maxPeriod = 0;
	numElements = 0;
	type = TT_EXPECTATION;
	numSamples = 0;

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
	sampleArray[numSamples] = TimeSample(time, state);
	numSamples++;
	my_adaboost.add_time(time, state > 0.5);
	return 0;
}

/*not required in incremental version*/
void CAdaboost::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	CFrelement fremen(0);
	fremen.init(maxPeriod, id, 1);
	for (int i = 0; i < numSamples; ++i) {
		fremen.add(sampleArray[i].t, sampleArray[i].v);
	}
	fremen.update(id);
	for (int i = 0; i < id; ++i) {
		my_adaboost.add_period(fremen.getPredictFrelements()[i].period);
	}
	my_adaboost.train();

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
	my_adaboost.print();
}

float CAdaboost::estimate(uint32_t time)
{
	return (my_adaboost.classify(time)+1)/2;
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
	my_adaboost.save(file, lossy);
	return 0;
}

int CAdaboost::load(FILE* file)
{
	my_adaboost.load(file);
	return 0;
}

int CAdaboost::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = type;
	my_adaboost.exportToArray(array, maxLen, pos);
	return pos;
}

int CAdaboost::importFromArray(double* array,int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_ADABOOST) std::cerr << "Error loading the model, type mismatch." << std::endl;
	my_adaboost.importFromArray(array, len, pos);
	return pos;
}

