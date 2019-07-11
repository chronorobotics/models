#ifndef CMOMENTS_H
#define CMOMENTS_H

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "CTemporal.h"

#define MAX_ID_LENGTH 100
	
using namespace std;

class CMoments: public CTemporal
{
	public:
		CMoments(int idd);
		~CMoments();

		void init(int iMaxPeriod, int elements, int numClasses);

		//adds a serie of measurements to the data
		int add(uint32_t time,float state);

		float estimate(uint32_t time);

		float predict(uint32_t time);

		void update(int maxOrder,unsigned int* times = NULL,float* signal = NULL,int length = 0);
		void print(bool verbose=true);

		int save(const char* name,bool lossy = false);
		int load(const char* name);

		int exportToArray(double* array,int maxLen);
		int importFromArray(double* array,int len);
		int save(FILE* file,bool lossy = false);
		int load(FILE* file);

	private:
		int id;
		float estimation;
		int positives, negatives;

		struct RightSide {
				double moment1_re;
				double moment1_im;
		};

		double pos_sum1_re;
		double pos_sum1_im;
		double neg_sum1_re;
		double neg_sum1_im;

		double pos_kappa;
		double pos_mu;
		double neg_kappa;
		double neg_mu;

		static int moment_f(const gsl_vector* x, void* params, gsl_vector* f);
};

#endif
