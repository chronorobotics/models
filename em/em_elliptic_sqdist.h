#ifndef EM_ELLIPTIC_SQDIST_H
#define EM_ELLIPTIC_SQDIST_H

#include <stdint.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "em_circular.h"

class EMEllipticSqdist : public EMCircular
{
	public:
		EMEllipticSqdist();
		EMEllipticSqdist(int cluster_count_, std::string filename);
		~EMEllipticSqdist();

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
		double maximisation();

		class Cluster {
			public:
				Cluster();
				Cluster(double xx_, double t0_, double aa_, double weight_);
				double xx;
				double t0;
				double aa;
				double weight;

				void print();
				void save(FILE* file, bool lossy = false);
				void load(FILE* file);
				void exportToArray(double* array, int maxLen, int& pos);
				void importFromArray(double* array, int len, int& pos);

				double density_at(double phase) const;
		};

		std::vector<Cluster> clusters;
		double timestamps_weight;
		std::ofstream test;
};

#endif // EM_ELLIPTIC_SQDIST_H
