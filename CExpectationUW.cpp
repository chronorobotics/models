#include <iostream>
#include <fstream>
#include <ctime>
#include "CExpectationUW.h"

CExpectationUW::CExpectationUW(int idd) :
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

void CExpectationUW::init(int iMaxPeriod, int elements, int numClasses)
{
	maxPeriod = iMaxPeriod;
	numElements = 1;
}

CExpectationUW::~CExpectationUW()
{
}

// adds new state observations at given times
int CExpectationUW::add(uint32_t time, float state)
{
	throw "We don't do this here";
}

int CExpectationUW::add_v(double x, double y, uint32_t time)
{
	model.add_value(x, y, time, 1);
	pedestrians++;
	return 0;
}

/*not required in incremental version*/
void CExpectationUW::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	model.train();
	std::cout << "Trained" << std::endl;
	print();
}

/*text representation of the fremen model*/
void CExpectationUW::print(bool verbose)
{
	model.print();
}

float CExpectationUW::estimate_v(double x, double y, uint32_t time)
{
	return model.get_density_at(x, y, time) * pedestrians;
}

float CExpectationUW::estimate(uint32_t time)
{
	throw "We don't do this here.";
}

float CExpectationUW::predict(uint32_t time)
{
	return estimate(time);
}

int CExpectationUW::save(const char* name, bool lossy)
{
	FILE* file = fopen(name,"w");
	save(file);
	fclose(file);
	return 0;
}

int CExpectationUW::load(const char* name)
{
	FILE* file = fopen(name,"r");
	load(file);
	fclose(file);
	return 0;
}


int CExpectationUW::save(FILE* file, bool lossy)
{
	fwrite(&pedestrians, sizeof(int), 1, file);
	model.save(file, lossy);
	return 0;
}

int CExpectationUW::load(FILE* file)
{
	fread(&pedestrians, sizeof(int), 1, file);
	model.load(file);
	return 0;
}

int CExpectationUW::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = TT_EXPECTATION;
	array[pos++] = pedestrians;
	model.exportToArray(array, maxLen, pos);
	return pos;
}

int CExpectationUW::importFromArray(double* array,int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_EXPECTATION) std::cerr << "Error loading the model, type mismatch." << std::endl;
	pedestrians = array[pos++];
	model.importFromArray(array, len, pos);
	return pos;
}

