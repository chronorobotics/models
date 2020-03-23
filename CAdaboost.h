#ifndef CADABOOST_H
#define CADABOOST_H

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include "CTemporal.h"

#include "adaboost/adaboost.h"
#include "em/em_gaussian.h"

class CAdaboost : public CTemporal
{
	public:
		CAdaboost(int idd);
		~CAdaboost();

		void init(int iMaxPeriod, int elements, int numClasses);

		//adds a serie of measurements to the data
		int add(uint32_t time,float state);


		//estimates the probability for the given times - using stored histogram
		float estimate(uint32_t time);

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
				TimeSample(long int t_, float v_);
				long int t;
				float v;
		};
		TimeSample sampleArray[1000000];
		int numSamples;
		long int min_time;
		long int max_time;

		std::vector<Adaboost> adaboosts;
		std::vector<double> means;
		EMGaussian class_model;
};

#endif // CADABOOST_H
