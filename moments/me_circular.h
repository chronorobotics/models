#ifndef ME_CIRCULAR_H
#define ME_CIRCULAR_H

#include "moment_estimator.h"

class MECircular : public MomentEstimator
{
	public:
		MECircular(int moment_count_);
		~MECircular();

		void add_point(uint32_t time);
		std::vector<double> get_moments() const;

		static double time_to_phase(uint32_t time);

		int get_moment_count() const;

	private:
		int moment_count;
};

#endif // ME_CIRCULAR_H
