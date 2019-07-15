#ifndef RIGHT_SIDE_H
#define RIGHT_SIDE_H

class MomentEstimator;

class RightSide {
	public:
		RightSide();
		RightSide(const MomentEstimator& me);
		~RightSide();
		double* moment_re;
		double* moment_im;
		int count;
};

#endif // RIGHT_SIDE_H
