#ifndef EM_CIRCULAR_H
#define EM_CIRCULAR_H

#include <vector>
#include <stdint.h>
#include <stdio.h>

class EMCircular
{
	public:
		EMCircular(int cluster_count_, double period_);

		virtual void train() = 0;
		virtual void add_time(uint32_t time) = 0;

		virtual void save(FILE* file, bool lossy = false) = 0;
		virtual void load(FILE* file) = 0;
		virtual void exportToArray(double* array, int maxLen, int& pos) = 0;
		virtual void importFromArray(double* array, int len, int& pos) = 0;

		virtual void print() = 0;

		virtual double get_density_at(uint32_t time) = 0;

	protected:
		int cluster_count;
		double period;
		static double time_to_phase(uint32_t time, double period_);

		class Timestamp {
			public:
				Timestamp(uint32_t time, int cluster_count, double period_);
				double phase;
				std::vector<double> alpha;
		};

		std::vector<Timestamp> timestamps;
};

#endif // EM_CIRCULAR_H
