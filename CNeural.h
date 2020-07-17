#ifndef CNEURAL_H
#define CNEURAL_H

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include "CTemporal.h"

#include "neural/neuronmodel.h"

class CNeural : public CTemporal
{
	public:
		CNeural(int idd);
		CNeural(int idd, int dimension_);
		~CNeural();

		void init(int iMaxPeriod, int elements, int numClasses);

		//adds a serie of measurements to the data
		int add(uint32_t time,float state);
		int add_v(uint32_t time, const std::vector<double>& state);

		//estimates the probability for the given times - using stored histogram
		float estimate(uint32_t time);
		std::vector<float> estimate_v(uint32_t time);

		//predicts the probability for the given times - using updated histogram
		float predict(uint32_t time);

		void update(int maxOrder, unsigned int* times = NULL, float* signal = NULL, int length = 0);
		void print(bool verbose=true);

		int save(const char* name, bool lossy = false);
		int load(const char* name);
		int save(FILE* file, bool lossy = false);
		int load(FILE* file);
		int exportToArray(double* array, int maxLen);
		int importFromArray(double* array, int len);

	private:
		int id;

		class TimeSample {
			public:
				TimeSample();
				TimeSample(long int t_, const std::vector<double>& v_);
				long int t;
				std::vector<double> v;
		};
		std::vector<TimeSample> sampleArray;
		int numSamples;
		int dimension;
		long int min_time;
		long int max_time;

		NeuronModel my_model;
};

#endif // CNEURAL_H
