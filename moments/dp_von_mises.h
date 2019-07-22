#ifndef DP_VON_MISES_H
#define DP_VON_MISES_H

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "density_params.h"

class DPVonMises : public DensityParams
{
	public:
		DPVonMises(int count_ = -1);
		~DPVonMises();
		std::vector<double> kappa;
		std::vector<double> mu;
		std::vector<double> weight;

		void calculate(RightSide& rs);
		double density_at(double phase);
		void print();
		void reset(int count_);

		int get_param_count() {
			return 3;
		}

		static double lnhyp(double x);
		static double hyp(double x);

		static int moment_f(const gsl_vector* x, void* params, gsl_vector* f);
};

#endif // DP_VON_MISES_H
