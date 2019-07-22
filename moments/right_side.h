#ifndef RIGHT_SIDE_H
#define RIGHT_SIDE_H

#include <vector>

class MomentEstimator;

class RightSide {
	public:
		RightSide();
		RightSide(const MomentEstimator& me);
		~RightSide();
		std::vector<double> moment_re;
		std::vector<double> moment_im;
		int count;
};

#endif // RIGHT_SIDE_H
