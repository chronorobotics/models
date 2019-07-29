#ifndef DP_SQDIST_H
#define DP_SQDIST_H

#include <vector>
#include <memory>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

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
				~EquationParams();
				std::vector<double> right_side;
				int cluster_count;
				int moment_count;
				gsl_integration_cquad_workspace* workplace;
		};

		class IntegrationParams {
			public:
				IntegrationParams(int u_, int v_, int n_);
				int u, v, n;
		};

		static int moment_f(const gsl_vector* x, void* params, gsl_vector* f);
		static double sin_integral_f(double x, void* params);
		static double cos_integral_f(double x, void* params);

		void importFromArray(double* array, int len, int& pos);
		int load(FILE* file);

		std::shared_ptr<MECircular> estimator;
};

#endif // DP_SQDIST_H
