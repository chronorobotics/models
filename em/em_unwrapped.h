#ifndef EM_UNWRAPPED_H
#define EM_UNWRAPPED_H

#include <vector>
#include "em_circular.h"

class EMUnwrapped : public EMCircular
{
	public:
	public:
		EMUnwrapped();
		EMUnwrapped(int cluster_count_);

		void train();
		void add_time(uint32_t time, double value);

		void save(FILE* file, bool lossy = false);
		void load(FILE* file);
		void exportToArray(double* array, int maxLen, int& pos);
		void importFromArray(double* array, int len, int& pos);

		void print();

		double get_density_at_d(double phase) const;

		/*std::vector<double> get_alpha_at(double value) const;
		std::vector<double> get_means() const;
		int get_cluster_count() const;*/

	private:
		void expectation();
		double maximisation();

		int cluster_count;

		class Cluster {
			public:
				Cluster();
				Cluster(double mu_, double sigma_, double weight_);
				double mu;
				double sigma;
				double weight;

				void print();
				void save(FILE* file, bool lossy = false);
				void load(FILE* file);
				void exportToArray(double* array, int maxLen, int& pos);
				void importFromArray(double* array, int len, int& pos);

				double density_at(double value) const;
		};

		std::vector<Cluster> clusters;
};

#endif // EM_UNWRAPPED_H
