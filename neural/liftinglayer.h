#ifndef LIFTINGLAYER_H
#define LIFTINGLAYER_H

#include "neuronlayer.h"

class LiftingLayer : public NeuronLayer
{
	public:
		LiftingLayer();
		LiftingLayer(std::vector<double> periods_, double t_min, double t_max);
		~LiftingLayer();

		virtual std::vector<double> forward(const std::vector<double>& x, int batch_size);
		virtual std::vector<double> backward(const std::vector<double>& dL_wrt_output);
		virtual void step(double rate);

		virtual void save(FILE* file, bool lossy = false);
		virtual void load(FILE* file);
		virtual void exportToArray(double* array, int maxLen, int& pos);
		virtual void importFromArray(double* array, int len, int& pos);

		virtual void print();

	private:
		std::vector<double> periods;
		double scale, offset;

		std::vector<double> last_x;
		int bs;
};

#endif // LIFTINGLAYER_H
