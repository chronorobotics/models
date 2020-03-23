#include <iostream>
#include <fstream>
#include "CMoments.h"

using namespace std;

CMoments::CMoments(int idd) :
	cluster_count(idd),
	estimators(),
	densities(),
	class_model(7)
{
	id=idd;
	measurements = 0;
	type = TT_MOMENTS;

	numSamples = 0;

	for (int i = class_model.get_cluster_count(); i; --i) {
		densities.push_back(std::move(DensityParams::create(this, DensityParams::SQDIST)));
		estimators.push_back(densities[densities.size()-1]->get_moment_estimator());
	}

	srandom(time(0));

	/*pos_density = DensityParams::create(this, DensityParams::SQDIST);
  neg_density = DensityParams::create(this, DensityParams::SQDIST);

	pos_estimator = pos_density->get_moment_estimator();
	neg_estimator = neg_density->get_moment_estimator();*/
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

int CMoments::get_cluster_count() const {
	return cluster_count;
}

/*double CMoments::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 86400.0f) / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}*/

// adds new state observations at given times
int CMoments::add(uint32_t time,float state)
{
	/*if (rand() % 100) {
    return 0;
  }
	sampleArray[numSamples] = TimeSample(time, state);
	numSamples++;

	if (state > 0.5) {
		pos_estimator->add_point(time, 1);
	} else {
		neg_estimator->add_point(time, 1);
	}
	measurements++;
	return 0;*/
	sampleArray[numSamples] = TimeSample(time, state);
	numSamples++;

	class_model.add_value(state);
	return 0;
}

void CMoments::update(int modelOrder, unsigned int* times, float* signal, int length)
{
	class_model.train();
	means = class_model.get_means();
	std::cout << numSamples << " " << means.size() << " " << densities.size() << std::endl;
	for (int i = 0; i < numSamples; ++i) {
		std::vector<double> alpha = class_model.get_alpha_at(sampleArray[i].v);
		for (int j = 0; j < alpha.size(); ++j) {
			estimators[j]->add_point(sampleArray[i].t, alpha[j]);
		}
	}

	for (int i = 0; i < densities.size(); ++i) {
		densities[i]->calculate();
	}

	//pos_density->calculate();
  //neg_density->calculate();

	ofstream myfile0("0.txt");
	ofstream myfile1("1.txt");
	for (int i = 0; i < numSamples; ++i) {
		myfile0 << sampleArray[i].t << " " << sampleArray[i].v << std::endl;
		myfile1 << sampleArray[i].t << " " << estimate(sampleArray[i].t) << std::endl;
	}
	myfile0.close();
	myfile1.close();
	print();
}

/*text representation of the fremen model*/
void CMoments::print(bool verbose)
{
	std::cout << "Model " << id << " Size: " << measurements << " " << std::endl;
	if (verbose) {
		for (int i = 0; i < means.size(); ++i) {
			std::cout << "mean " << means[i] << ":";
			densities[i]->print();
			std::cout << std::endl;
		}
	}
}

float CMoments::estimate(uint32_t time)
{
	/*float pd = pos_density->density_at(time);
	//float nd = neg_density->density_at(time);

	//return pd / (pd + nd);
	return pd;*/
	double s = 0;
	double t = 0;
	for (int i = 0; i < means.size(); ++i) {
		double density = densities[i]->density_at(time);
		s += means[i] * density;
		t += density;
	}

	return s/t;
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
/*	pos_density->save(file, lossy);
	neg_density->save(file, lossy);*/
	int size = densities.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < densities.size(); ++i) {
		fwrite(&(means[i]), sizeof(double), 1, file);
		densities[i]->save(file, lossy);
	}
	return 0;
}

int CMoments::load(FILE* file)
{
/*	pos_density = DensityParams::load(this, file);
	neg_density = DensityParams::load(this, file);*/
	int size;
	fread(&size, sizeof(int), 1, file);
	means.resize(size);
	densities.resize(size);
	for (int i = 0; i < densities.size(); ++i) {
		fread(&(means[i]), sizeof(double), 1, file);
		densities[i] = DensityParams::load(this, file);
	}
	return 0;
}


int CMoments::exportToArray(double* array, int maxLen)
{
	/*int pos = 0;
	array[pos++] = type;
	pos_density->exportToArray(array, maxLen, pos);
	neg_density->exportToArray(array, maxLen, pos);*/
	int pos = 0;
	array[pos++] = type;
	array[pos++] = densities.size();
	for (int i = 0; i < densities.size(); ++i) {
		array[pos++] = means[i];
		densities[i]->exportToArray(array, maxLen, pos);
	}
	return pos;
}

int CMoments::importFromArray(double* array, int len)
{
	/*int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_NONE) std::cerr << "Error loading the model, type mismatch." << std::endl;
	pos_density = DensityParams::importFromArray(this, array, len, pos);
	neg_density = DensityParams::importFromArray(this, array, len, pos);*/
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_MOMENTS) std::cerr << "Error loading the model, type mismatch." << std::endl;
	int size = array[pos++];
	means.resize(size);
	densities.resize(size);
	for (int i = 0; i < densities.size(); ++i) {
		means[i] = array[pos++];
		densities[i] = DensityParams::importFromArray(this, array, len, pos);
	}
	return pos;
}
