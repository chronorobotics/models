#ifndef MINMAXLAYER_H
#define MINMAXLAYER_H

#include "neuronlayer.h"

class MinMaxLayer : public NeuronLayer
{
	public:
		MinMaxLayer();
		MinMaxLayer(int size_);
		~MinMaxLayer();

		virtual std::vector<double> forward(const std::vector<double>& x, int batch_size);
		virtual std::vector<double> backward(const std::vector<double>& dL_wrt_output);
		virtual void step(double rate);

		virtual void save(FILE* file, bool lossy = false);
		virtual void load(FILE* file);
		virtual void exportToArray(double* array, int maxLen, int& pos);
		virtual void importFromArray(double* array, int len, int& pos);

		virtual void print();

	private:
		int size;
		std::vector<double> last_x;
		int bs;
};

#endif // MINMAXLAYER_H
