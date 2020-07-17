#ifndef EM_GAUSS_VONMISES_H
#define EM_GAUSS_VONMISES_H

#include <stdio.h>
#include <vector>
#include <stdint.h>

class EMGaussVonMises
{
	public:
		EMGaussVonMises();
		EMGaussVonMises(int cluster_count_);

		void train();
		void add_value(double x, double y, uint32_t time, double weight);

		void save(FILE* file, bool lossy = false);
		void load(FILE* file);
		void exportToArray(double* array, int maxLen, int& pos);
		void importFromArray(double* array, int len, int& pos);

		void print();

		double get_density_at(double x, double y, uint32_t time) const;
		int get_cluster_count() const;

	private:
		void expectation();
		double maximisation();

		static double time_to_phase(uint32_t time);

		int cluster_count;
		double total_weight;

		class Point {
			public:
				Point(double x_, double y_, uint32_t time, int cluster_count, double weight_);
				double x, y, phase;
				std::vector<double> alpha;
				double weight;
		};

		class Cluster {
			public:
				Cluster();
				Cluster(double ex_, double ey_, double cxx_, double cyy_, double cxy_, double det_, double tmu_, double tkappa_, double weight_);
				double ex, ey, cxx, cyy, cxy, det, tmu, tkappa;
				double weight;

				void print();
				void save(FILE* file, bool lossy = false);
				void load(FILE* file);
				void exportToArray(double* array, int maxLen, int& pos);
				void importFromArray(double* array, int len, int& pos);

				double density_at(double x, double y, double phase) const;
				void estimate_from_mean(double re, double im, bool keep_kappa);

				static double mean_f (double x, void* params);
				static double mean_df (double x, void* params);
				static void mean_fdf (double x, void* params, double* y, double* dy);
		};

		std::vector<Cluster> clusters;
		std::vector<Point> points;
};

#endif // EM_GAUSS_VONMISES_H
