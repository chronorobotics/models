#ifndef MOMENT_ESTIMATOR_H
#define MOMENT_ESTIMATOR_H

#include <stdint.h>
#include <vector>

class MomentEstimator {
	public:
		MomentEstimator();
		~MomentEstimator();

		virtual void add_point(uint32_t time, double weight) = 0;
		virtual std::vector<double> get_moments() const = 0;

	protected:
		std::vector<double> data;
		double count;
};

#endif // MOMENT_ESTIMATOR_H
