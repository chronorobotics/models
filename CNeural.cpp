#include <iostream>
#include <fstream>
#include <ctime>
#include <limits>
#include "CNeural.h"
#include "CFrelement.h"

#include "neural/liftinglayer.h"
#include "neural/linearlayer.h"
#include "neural/logreglayer.h"
#include "neural/minmaxlayer.h"

CNeural::CNeural(int idd) :
	min_time(std::numeric_limits<long int>::max()),
	max_time(std::numeric_limits<long int>::min()),
	my_model(1)
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

void CNeural::init(int iMaxPeriod, int elements, int numClasses)
{
	maxPeriod = iMaxPeriod;
	numElements = 1;
}

CNeural::~CNeural()
{
}


CNeural::TimeSample::TimeSample() :
	t(),
	v()
{
}

CNeural::TimeSample::TimeSample(long t_, float v_) :
	t(t_),
	v(v_)
{
}

// adds new state observations at given times
int CNeural::add(uint32_t time, float state)
{
	sampleArray[numSamples] = TimeSample(time, state);
	numSamples++;
	if (time < min_time) {
		min_time = time;
	}
	if (time > max_time) {
		max_time = time;
	}
	return 0;
}

/*not required in incremental version*/
void CNeural::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	CFrelement fremen(0);
	fremen.init(maxPeriod, id, 1);
	for (int i = 0; i < numSamples; ++i) {
		fremen.add(sampleArray[i].t, sampleArray[i].v);
	}
	fremen.update(id);
	std::vector<double> periods;
	for (int i = 0; i < id; ++i) {
		periods.push_back(fremen.getPredictFrelements()[i].period);
	}
	int id1 = id;
	int id2 = id1*id1;
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LiftingLayer(periods, min_time, max_time)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LogregLayer(id1, 7)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new MinMaxLayer(7*id1)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LinearLayer(49*id2, 3*id1)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new MinMaxLayer(3*id1)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LinearLayer(9*id2, 5)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new MinMaxLayer(5)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LinearLayer(25, 1)));
	//my_model.add_layer(std::shared_ptr<NeuronLayer>(new LinearLayer(7*id1, 1)));

	double err;
	int small = 0;
	bool smallstep = false;
	ofstream myfileerr("err.txt");
	for (int i = 0; i < 1000000; ++i) {
		std::cout << "iteration " << i << std::endl;
		int batch_size = 100;
		std::vector<double> x;
		std::vector<double> y;
		x.reserve(batch_size);
		y.reserve(batch_size);
		for (int j = 0; j < batch_size; ++j) {
			int r = rand() % numSamples;
			x.push_back(sampleArray[r].t);
			y.push_back(sampleArray[r].v);
		}
		err = my_model.forward(x, y);
		if (err < 0.12) {
			small++;
			if (small > 30) {
				smallstep = true;
			}
		} else {
			small = 0;
		}
		myfileerr << err << std::endl;
		my_model.backward();
		my_model.step(smallstep ? 1E-3 : 1E-2);
		/*if (i == 50000) {
			my_model.simplify();
		}*/
	}
	myfileerr.close();

	print();

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
void CNeural::print(bool verbose)
{
	my_model.print();
}

float CNeural::estimate(uint32_t time)
{
	std::vector<double> x;
	x.push_back(time);
	std::vector<double> y = my_model.forward(x);
	return std::max(y[0], 0.0);
}

float CNeural::predict(uint32_t time)
{
	return estimate(time);
}

int CNeural::save(const char* name, bool lossy)
{
	FILE* file = fopen(name,"w");
	save(file);
	fclose(file);
	return 0;
}

int CNeural::load(const char* name)
{
	FILE* file = fopen(name,"r");
	load(file);
	fclose(file);
	return 0;
}


int CNeural::save(FILE* file, bool lossy)
{
	my_model.save(file, lossy);
	return 0;
}

int CNeural::load(FILE* file)
{
	my_model.load(file);
	return 0;
}

int CNeural::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = type;
	my_model.exportToArray(array, maxLen, pos);
	return pos;
}

int CNeural::importFromArray(double* array,int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_NEURAL) std::cerr << "Error loading the model, type mismatch." << std::endl;
	my_model.importFromArray(array, len, pos);
	return pos;
}

