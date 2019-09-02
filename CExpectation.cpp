#include <iostream>
#include <fstream>
#include <ctime>
#include "CExpectation.h"

CExpectation::CExpectation(int idd) :
	pedestrians(0),
	model(idd)
{
	id=idd;
	firstTime = -1;
	lastTime = -1;
	measurements = 0;
	maxPeriod = 0;
	numElements = 0;
	type = TT_EXPECTATION;

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

// adds new state observations at given times
int CExpectation::add(uint32_t time, float state)
{
	throw "We don't do this here";
}

int CExpectation::add_v(double x, double y, uint32_t time)
{
	model.add_value(x, y, time, 1);
	pedestrians++;
	return 0;
}

/*not required in incremental version*/
void CExpectation::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	model.train();
	std::cout << "Trained" << std::endl;
}

/*text representation of the fremen model*/
void CExpectation::print(bool verbose)
{
	model.print();
}

float CExpectation::estimate_v(double x, double y, uint32_t time)
{
	return model.get_density_at(x, y, time);
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
	fwrite(&pedestrians, sizeof(int), 1, file);
	model.save(file, lossy);
	return 0;
}

int CExpectation::load(FILE* file)
{
	fread(&pedestrians, sizeof(int), 1, file);
	model.load(file);
	return 0;
}

int CExpectation::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = TT_EXPECTATION;
	array[pos++] = pedestrians;
	model.exportToArray(array, maxLen, pos);
	return pos;
}

int CExpectation::importFromArray(double* array,int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_EXPECTATION) std::cerr << "Error loading the model, type mismatch." << std::endl;
	pedestrians = array[pos++];
	model.importFromArray(array, len, pos);
	return pos;
}

