#ifndef NEURONLAYER_H
#define NEURONLAYER_H

#include <stdio.h>
#include <iostream>
#include <vector>

class NeuronLayer
{
	public:
		NeuronLayer();
		~NeuronLayer();

		enum NeuronLayerType {
			NL_LIFTING,
			NL_LOGREG,
			NL_LINEAR,
			NL_MINMAX
		};

		virtual std::vector<double> forward(const std::vector<double>& x, int batch_size) = 0;
		virtual std::vector<double> backward(const std::vector<double>& dL_wrt_output) = 0;
		virtual void step(double rate) = 0;
		virtual void simplify();

		virtual void save(FILE* file, bool lossy = false) = 0;
		virtual void load(FILE* file) = 0;
		virtual void exportToArray(double* array, int maxLen, int& pos) = 0;
		virtual void importFromArray(double* array, int len, int& pos) = 0;

		virtual void print() = 0;

	protected:
		void random_init(std::vector<double>& vec, double scale);
};

#endif // NEURONLAYER_H
