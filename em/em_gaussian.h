#ifndef EM_GAUSSIAN_H
#define EM_GAUSSIAN_H

#include <vector>

class EMGaussian
{
	public:
	public:
		EMGaussian(int cluster_count_);

		void train();
		void add_value(double value);

		/*void save(FILE* file, bool lossy = false);
		void load(FILE* file);
		void exportToArray(double* array, int maxLen, int& pos);
		void importFromArray(double* array, int len, int& pos);*/

		void print();

		std::vector<double> get_alpha_at(double value) const;
		std::vector<double> get_means() const;
		int get_cluster_count() const;

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
				/*void save(FILE* file, bool lossy = false);
				void load(FILE* file);
				void exportToArray(double* array, int maxLen, int& pos);
				void importFromArray(double* array, int len, int& pos);*/

				double density_at(double value) const;
		};

		class Measurement {
			public:
				Measurement(double value_, int cluster_count);
				double value;
				std::vector<double> alpha;
		};

		std::vector<Cluster> clusters;
		std::vector<Measurement> measurements;
};

#endif // EM_GAUSSIAN_H
