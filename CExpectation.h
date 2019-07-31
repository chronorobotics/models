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

class CExpectation : public CTemporal
{
	public:
		CExpectation(int idd);
		~CExpectation();

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

		static const int cluster_count = 1;
		static double time_to_phase(uint32_t time);

		void expectation();
		void maximisation();

		class Cluster {
			public:
				Cluster();
				Cluster(double kappa_, double mu_, double weight_);
				double kappa;
				double mu;
				double weight;

				void print();
				void save(FILE* file, bool lossy = false);
				void load(FILE* file);
				void exportToArray(double* array, int maxLen, int& pos);
				void importFromArray(double* array, int len, int& pos);

				double density_at(double phase) const;
				void estimate_from_mean(double re, double im);

				static double mean_f (double x, void* params);
				static double mean_df (double x, void* params);
				static void mean_fdf (double x, void* params, double* y, double* dy);
		};

		class Timestamp {
			public:
				Timestamp(uint32_t time, int cluster_count);
				double phase;
				std::vector<double> alpha;
		};

		std::vector<Cluster> clusters;
		std::vector<Timestamp> timestamps;

		class TimeSample {
			public:
				TimeSample();
				TimeSample(long int t_, float v_);
				long int t;
				float v;
		};
		TimeSample sampleArray[1000000];
		int numSamples;
};

#endif // CEXPECTATION_H
