#ifndef CMOMENTS_H
#define CMOMENTS_H

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <memory>
#include "CTemporal.h"

#include "moments/moment_estimator.h"
#include "moments/density_params.h"

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

		int get_cluster_count() const;

	private:
		int id;
		float estimation;
		int cluster_count;

		class TimeSample {
			public:
				TimeSample();
				TimeSample(long int t_, float v_);
				long int t;
				float v;
		};

		MomentEstimator* pos_estimator;
		MomentEstimator* neg_estimator;
		std::unique_ptr<DensityParams> pos_density;
		std::unique_ptr<DensityParams> neg_density;

		TimeSample sampleArray[1000000];
		int numSamples;
		//double time_to_phase(uint32_t time);
};

#endif
