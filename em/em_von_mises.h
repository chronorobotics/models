#ifndef EM_VON_MISES_H
#define EM_VON_MISES_H

#include <stdint.h>
#include <vector>

#include "em_circular.h"

class EMVonMises : public EMCircular
{
	public:
		EMVonMises();
		EMVonMises(int cluster_count_);

		void train();
		void add_time(uint32_t time, double value);

		void save(FILE* file, bool lossy = false);
		void load(FILE* file);
		void exportToArray(double* array, int maxLen, int& pos);
		void importFromArray(double* array, int len, int& pos);

		void print();

		double get_density_at(uint32_t time);

	private:
		void expectation();
		double maximisation(bool keep_kappa);

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
				void estimate_from_mean(double re, double im, bool keep_kappa);

				static double mean_f (double x, void* params);
				static double mean_df (double x, void* params);
				static void mean_fdf (double x, void* params, double* y, double* dy);
		};

		std::vector<Cluster> clusters;
		double timestamps_weight;
};

#endif // EM_VON_MISES_H
