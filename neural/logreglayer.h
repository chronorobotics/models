#ifndef LOGREGLAYER_H
#define LOGREGLAYER_H

#include "neuronlayer.h"

class LogregLayer : public NeuronLayer
{
	public:
		LogregLayer();
		LogregLayer(int channels_, int sets_);
		~LogregLayer();

		virtual std::vector<double> forward(const std::vector<double>& x, int batch_size);
		virtual std::vector<double> backward(const std::vector<double>& dL_wrt_output);
		virtual void step(double rate);

		virtual void save(FILE* file, bool lossy = false);
		virtual void load(FILE* file);
		virtual void exportToArray(double* array, int maxLen, int& pos);
		virtual void importFromArray(double* array, int len, int& pos);

		virtual void print();

	private:
		int channels;
		int sets;
		std::vector<double> xw;
		std::vector<double> yw;
		std::vector<double> zw;

		std::vector<double> last_x;
		int bs;
		std::vector<double> dL_wrt_xw;
		std::vector<double> dL_wrt_yw;
		std::vector<double> dL_wrt_zw;
};

#endif // LOGREGLAYER_H
