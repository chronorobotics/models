#ifndef NEURONMODEL_H
#define NEURONMODEL_H

#include <memory>
#include <vector>

#include "squareerror.h"

class NeuronLayer;

class NeuronModel
{
	public:
		NeuronModel(int output_size_);

		std::vector<double> forward(const std::vector<double>& x_);
		double forward(const std::vector<double>& x_, const std::vector<double>& y);
		void backward();
		void step(float rate);
		void simplify();

		void add_layer(std::shared_ptr<NeuronLayer> layer);

		void save(FILE* file, bool lossy = false);
		void load(FILE* file);
		void exportToArray(double* array, int maxLen, int& pos);
		void importFromArray(double* array, int len, int& pos);

		void print();

	private:
		std::vector<std::shared_ptr<NeuronLayer> > layers;
		SquareError se;
		int output_size;
};

#endif // NEURONMODEL_H
