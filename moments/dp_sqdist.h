#ifndef DP_SQDIST_H
#define DP_SQDIST_H

#include <vector>
#include <memory>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "density_params.h"

class MECircular;

class DPSqdist : public DensityParams
{
	public:
	public:
		DPSqdist(CMoments* parent_, int count_ = -1);
		~DPSqdist();
		std::vector<double> xx;
		std::vector<double> yy;
		std::vector<double> weight;

		virtual MomentEstimator* get_moment_estimator();
		void calculate();
		double density_at(uint32_t time);
		void print();
		void reset(int count_);

		void exportToArray(double* array, int maxLen, int& pos);
		int save(FILE* file, bool lossy = false);

	protected:
		class EquationParams {
			public:
				EquationParams(std::vector<double> right_side_, int cluster_count_, int moment_count_);
				std::vector<double> right_side;
				int cluster_count;
				int moment_count;
		};

		static int moment_f(const gsl_vector* x, void* params, gsl_vector* f);
		static void power(double& re, double& im, int n);

		void importFromArray(double* array, int len, int& pos);
		int load(FILE* file);

		std::shared_ptr<MECircular> estimator;
};

#endif // DP_SQDIST_H
