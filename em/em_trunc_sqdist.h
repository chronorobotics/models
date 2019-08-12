#ifndef EM_TRUNC_SQDIST_H
#define EM_TRUNC_SQDIST_H

#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

class EMTruncSqdist
{
	public:
		EMTruncSqdist(int cluster_count_);

		void train();
		void add_time(uint32_t time);

		void save(FILE* file, bool lossy = false);
		void load(FILE* file);
		void exportToArray(double* array, int maxLen, int& pos);
		void importFromArray(double* array, int len, int& pos);

		void print();

		double get_density_at(uint32_t time);

	private:
		void expectation();
		double maximisation();

		int cluster_count;
		static double time_to_phase(uint32_t time);

		class Cluster {
			public:
				Cluster();
				Cluster(double xx_, double yy_, double theta_, double weight_);
				double xx;
				double yy;
				double theta;
				double weight;

				double corrective;
				double trunc;

				void print();
				void save(FILE* file, bool lossy = false);
				void load(FILE* file);
				void exportToArray(double* array, int maxLen, int& pos);
				void importFromArray(double* array, int len, int& pos);

				void precalc();
				void estimate(double m1r, double m1i, double m2r, double m2i);
				double density_at(double phase) const;

			private:
				class RightSide {
					public:
						RightSide(double m1r, double m1i, double m2r, double m2i);
						double m1, m2;
				};

				static int estimation_f(const gsl_vector* X, void* params, gsl_vector* f);
		};

		class Timestamp {
			public:
				Timestamp(uint32_t time, int cluster_count);
				double phase;
				std::vector<double> alpha;
		};

		std::vector<Cluster> clusters;
		std::vector<Timestamp> timestamps;
};

#endif // EM_TRUNC_SQDIST_H
