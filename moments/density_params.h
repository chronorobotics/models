#ifndef DENSITY_PARAMS_H
#define DENSITY_PARAMS_H

class CMoments;
class RightSide;

class DensityParams {
	public:
		DensityParams(CMoments* parent_, int count_ = -1);
		~DensityParams();
		int count;

		virtual int get_param_count() = 0;

		virtual void calculate(RightSide& rs) = 0;
		virtual double density_at(double phase) = 0;
		virtual void print() = 0;
		virtual void reset(int count_) = 0;

	protected:
		CMoments* parent;
};

#endif // DENSITY_PARAMS_H
