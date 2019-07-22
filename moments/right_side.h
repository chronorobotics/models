#ifndef RIGHT_SIDE_H
#define RIGHT_SIDE_H

#include <vector>

class CMoments;
class MomentEstimator;

class RightSide {
	public:
		RightSide(CMoments* parent_);
		RightSide(const MomentEstimator& me);
		~RightSide();
		std::vector<double> moment_re;
		std::vector<double> moment_im;
		int count;

	private:
		CMoments* parent;
};

#endif // RIGHT_SIDE_H
