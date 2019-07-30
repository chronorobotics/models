#ifndef DP_UNIFORM_H
#define DP_UNIFORM_H

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "density_params.h"

class MECircular;

class DPUniform : public DensityParams
{
	public:
		DPUniform(CMoments* parent_, int count_ = -1);
		~DPUniform();
		std::vector<double> start;
		std::vector<double> end;
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

		static double normalise_angle(double x);
		static double arc_length(double a, double b);

		static int moment_f(const gsl_vector* x, void* params, gsl_vector* f);

		void importFromArray(double* array, int len, int& pos);
		int load(FILE* file);

		std::shared_ptr<MECircular> estimator;
};

#endif // DP_UNIFORM_H
