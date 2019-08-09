#include <iostream>
#include <fstream>
#include <ctime>
#include "CExpectation.h"

CExpectation::CExpectation(int idd) :
	positive(idd),
	negative(idd)
{
	id=idd;
	firstTime = -1;
	lastTime = -1;
	measurements = 0;
	maxPeriod = 0;
	numElements = 0;
	type = TT_EXPECTATION;
	numSamples = 0;

	positives = 0;
	negatives = 0;

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

	if (state > 0.5) {
		positive.add_time(time);
		positives++;
	} else {
		negative.add_time(time);
		negatives++;
	}
	return 0;
}

/*not required in incremental version*/
void CExpectation::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	positive.train();
	negative.train();

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
	std::cout << "Positive:";
	positive.print();
	std::cout << std::endl << "Negative:";
	negative.print();
	std::cout << std::endl;
}

float CExpectation::estimate(uint32_t time)
{
	double pd = positive.get_density_at(time) * positives;
	double nd = negative.get_density_at(time) * negatives;

	return pd / (pd + nd);
	return pd;
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
	positive.save(file, lossy);
	negative.save(file, lossy);
	fwrite(&positives, sizeof(int), 1, file);
	fwrite(&negatives, sizeof(int), 1, file);
	return 0;
}

int CExpectation::load(FILE* file)
{
	positive.load(file);
	negative.load(file);
	fread(&positives, sizeof(int), 1, file);
	fread(&negatives, sizeof(int), 1, file);
	return 0;
}

int CExpectation::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = type;
	positive.exportToArray(array, maxLen, pos);
	negative.exportToArray(array, maxLen, pos);
	array[pos++] = positives;
	array[pos++] = negatives;
	return pos;
}

int CExpectation::importFromArray(double* array,int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_MEAN) std::cerr << "Error loading the model, type mismatch." << std::endl;
	positive.importFromArray(array, len, pos);
	negative.importFromArray(array, len, pos);
	positives = array[pos++];
	negatives = array[pos++];
	update(0);
	return pos;
}

