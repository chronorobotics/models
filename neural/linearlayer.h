#ifndef LINEARLAYER_H
#define LINEARLAYER_H

#include "neuronlayer.h"

class LinearLayer : public NeuronLayer
{
	public:
		LinearLayer();
		LinearLayer(int in, int out);
		~LinearLayer();

		virtual std::vector<double> forward(const std::vector<double>& x, int batch_size);
		virtual std::vector<double> backward(const std::vector<double>& dL_wrt_output);
		virtual void step(double rate);
		virtual void simplify();

		virtual void save(FILE* file, bool lossy = false);
		virtual void load(FILE* file);
		virtual void exportToArray(double* array, int maxLen, int& pos);
		virtual void importFromArray(double* array, int len, int& pos);

		virtual void print();

	private:
		int rows;
		int cols;
		std::vector<double> W;

		std::vector<double> last_x;
		int bs;
		std::vector<double> dL_wrt_W;
};

#endif // LINEARLAYER_H
