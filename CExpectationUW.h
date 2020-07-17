#ifndef CExpectationUW_H
#define CExpectationUW_H

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include "CTemporal.h"

#include "em/em_sqdist.h"
#include "em/em_von_mises.h"
#include "em/em_gaussian.h"
#include "em/em_unwrapped.h"

class CExpectationUW : public CTemporal
{
	public:
		CExpectationUW(int idd);
		~CExpectationUW();

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

		/*EMEllipticSqdist positive;
		EMEllipticSqdist negative;*/
		std::vector<EMUnwrapped> models;
		std::vector<double> means;
		EMGaussian class_model;
};

#endif // CExpectationUW_H
