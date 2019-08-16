#ifndef CEXPECTATION_H
#define CEXPECTATION_H

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include "CTemporal.h"

#include "em/em_sqdist.h"
#include "em/em_gaussian.h"
#include "em/em_multi_alternative.h"

class CExpectation : public CTemporal
{
	public:
		CExpectation(int idd);
		CExpectation(int idd, int dimension_);
		~CExpectation();

		void init(int iMaxPeriod, int elements, int numClasses);

		//adds a serie of measurements to the data
		int add(uint32_t time, float state);
		int add_v(uint32_t time, std::vector<bool> state);

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
				TimeSample(long int t_, std::vector<bool> v_);
				long int t;
				std::vector<bool> v;
		};
		TimeSample sampleArray[1000000];
		int numSamples;
		int dimension;

		/*EMSqdist positive;
		EMSqdist negative;*/
		std::vector<EMSqdist> models;
		std::vector<std::vector<double> > means;
		EMMultiAlternative class_model;
};

#endif // CEXPECTATION_H
