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

#include "em/em_gauss_vonmises.h"
#include "em/em_gauss_sqdist.h"
#include "em/em_gauss_unwrapped.h"

class CExpectationUW : public CTemporal
{
	public:
		CExpectationUW(int idd);
		~CExpectationUW();

		void init(int iMaxPeriod, int elements, int numClasses);

		//adds a serie of measurements to the data
		int add(uint32_t time, float state);
		int add_v(double x, double y, uint32_t time);

		//estimates the probability for the given times - using stored histogram
		float estimate(uint32_t time);
		float estimate_v(double x, double y, uint32_t time);

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

		int pedestrians;
		EMGaussUnwrapped model;
};

#endif // CExpectationUW_H
