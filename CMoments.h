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
		static const int moment_count = 12;

		class MomentEstimator {
			public:
				MomentEstimator();
				~MomentEstimator();
				double* sum_re;
				double* sum_im;
				int count;

				void add_point(double phase);
		};

		class RightSide {
			public:
				RightSide();
				RightSide(const MomentEstimator& me);
				~RightSide();
				double* moment_re;
				double* moment_im;
				int count;
		};

		class TimeSample {
			public:
				TimeSample();
				TimeSample(long int t_, float v_);
				long int t;
				float v;
		};

		class DensityParams {
			public:
				DensityParams(int count_ = -1);
				~DensityParams();
				double* kappa;
				double* mu;
				double* weight;
				int count;

				void calculate(RightSide& rs);
				double density_at(double phase);
				void print();
				void reset(int count_);
		};

		MomentEstimator pos_estimator;
		MomentEstimator neg_estimator;
		DensityParams pos_density;
		DensityParams neg_density;

		TimeSample sampleArray[1000000];
		int numSamples;
		double time_to_phase(uint32_t time);

		static int moment_f(const gsl_vector* x, void* params, gsl_vector* f);
		static int moment_df(const gsl_vector* x, void* params, gsl_matrix* J);
		static int moment_fdf (const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);
};

#endif
