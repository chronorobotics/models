#ifndef DP_DOUBLE_VON_MISES_H
#define DP_DOUBLE_VON_MISES_H

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "density_params.h"
#include "me_bicircular.h"

class DPDoubleVonMises : public DensityParams
{
	public:
		DPDoubleVonMises(CMoments* parent_, int count_ = -1);
		~DPDoubleVonMises();
		std::vector<double> kappa;
		std::vector<double> mu;
		std::vector<double> lambda;
		std::vector<double> nu;
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
				EquationParams(std::vector<double> right_side_, const std::vector<MEBicircular::Index>& indices_,
											 int cluster_count_, int moment_count_);
				std::vector<double> right_side;
				const std::vector<MEBicircular::Index>& indices;
				int cluster_count;
				int moment_count;
				int min_kappa;
				int min_lambda;
		};

		static double lnhyp(double x, double min_kappa);
		static double hyp(double x);

		static int moment_f(const gsl_vector* x, void* params, gsl_vector* f);

		void importFromArray(double* array, int len, int& pos);
		int load(FILE* file);

		std::shared_ptr<MEBicircular> estimator;
};

#endif // DP_DOUBLE_VON_MISES_H
