#ifndef ME_BICIRCULAR_H
#define ME_BICIRCULAR_H

#include "moment_estimator.h"

class MEBicircular : public MomentEstimator
{
	public:
		class Phase {
			public:
				Phase();
				Phase(uint32_t time);
				Phase(double f1_, double f2_);

				double f1;
				double f2;
		};

		class Index {
			public:
				Index();
				Index(int i1_, int i2_);

				int i1;
				int i2;
		};

		MEBicircular(int moment_count_);

		void add_point(uint32_t time);
		std::vector<double> get_moments() const;

		const std::vector<Index>& get_moment_indices() const;

	private:
		int moment_count;
		std::vector<Index> moment_indices;
};

#endif // ME_BICIRCULAR_H
