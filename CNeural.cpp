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
	throw "new WeDontDoThisHereException()";
}

CNeural::CNeural(int idd, int dimension_) :
	min_time(std::numeric_limits<long int>::max()),
	max_time(std::numeric_limits<long int>::min()),
	my_model(dimension_)
{
	id=idd;
	firstTime = -1;
	lastTime = -1;
	measurements = 0;
	maxPeriod = 0;
	numElements = 0;
	type = TT_EXPECTATION;
	numSamples = 0;
	dimension = dimension_;

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

CNeural::TimeSample::TimeSample(long t_, const std::vector<double>& v_) :
	t(t_),
	v(v_)
{
}

// adds new state observations at given times
int CNeural::add(uint32_t time, float state)
{
	throw "new WeDontDoThisHereException()";
}

int CNeural::add_v(uint32_t time, const std::vector<double>& state)
{
	sampleArray.push_back(TimeSample(time, state));
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
	/*CFrelement fremen(0);
	fremen.init(maxPeriod, id, 1);
	for (int i = 0; i < numSamples; ++i) {
		double v = 0;
		for (int j = 0; j < sampleArray[i].v.size(); ++j) {
			v += sampleArray[i].v[j] * tan(j);
		}
		fremen.add(sampleArray[i].t, v);
	}
	fremen.update(id);
	std::vector<double> periods;
	for (int i = 0; i < id; ++i) {
		periods.push_back(fremen.getPredictFrelements()[i].period);
	}
	int id1 = id;
	int id2 = id1*id1;
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LiftingLayer(periods, min_time, max_time)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LogregLayer(id1, 3)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new MinMaxLayer(3*id1)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LinearLayer(9*id2, id1)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new MinMaxLayer(id1)));
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LinearLayer(id2, dimension)));*/
	//my_model.add_layer(std::shared_ptr<NeuronLayer>(new LinearLayer(7*id1, 1)));

	FILE* f = fopen("neural_model", "r");
	my_model.load(f);
	fclose(f);
	my_model.pop_layer();
	my_model.add_layer(std::shared_ptr<NeuronLayer>(new LinearLayer(25, dimension)));
	my_model.set_output_size(dimension);

	double err;
	ofstream myfileerr("err.txt");
	for (int i = 0; i < 200; ++i) {
		std::cerr << "iteration " << i << std::endl;
		int batch_size = 30;
		std::vector<double> x;
		std::vector<double> y;
		x.reserve(batch_size);
		y.reserve(batch_size);
		for (int j = 0; j < batch_size; ++j) {
			int r = rand() % numSamples;
			x.push_back(sampleArray[r].t);
			y.insert(y.end(), sampleArray[r].v.begin(), sampleArray[r].v.end());
		}
		err = my_model.forward(x, y);
		myfileerr << err << std::endl;
		my_model.backward1();
		my_model.step1(1E-1);
		/*if (i == 50000) {
			my_model.simplify();
		}*/
	}
	myfileerr.close();

	print();
}

/*text representation of the fremen model*/
void CNeural::print(bool verbose)
{
	my_model.print();
}

float CNeural::estimate(uint32_t time)
{
	throw "new WeDontDoThisHereException()";
}

std::vector<float> CNeural::estimate_v(uint32_t time)
{
	std::vector<double> x;
	x.push_back(time);
	std::vector<double> y = my_model.forward(x);
	std::vector<float> result;
	for (int i = 0; i < y.size(); ++i) {
		result.push_back(y[i]);
	}
	return result;
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

