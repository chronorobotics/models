#ifndef DENSITY_PARAMS_H
#define DENSITY_PARAMS_H

#include <memory>

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

		virtual void exportToArray(double* array, int maxLen, int& pos) = 0;
		virtual int save(FILE* file, bool lossy = false) = 0;

		enum Distribution {
			VON_MISES = 0
		};

		static std::unique_ptr<DensityParams> create(CMoments* parent, Distribution dist, int count = -1);
		static std::unique_ptr<DensityParams> importFromArray(CMoments* parent, double* array, int len, int& pos);
		static std::unique_ptr<DensityParams> load(CMoments* parent, FILE* file);

	protected:
		CMoments* parent;

		virtual void importFromArray(double* array, int len, int& pos) = 0;
		virtual int load(FILE* file) = 0;
};

#endif // DENSITY_PARAMS_H
