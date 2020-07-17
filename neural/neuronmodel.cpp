#include <stdio.h>
#include <iostream>

#include "neuronlayer.h"
#include "neuronmodel.h"

#include "liftinglayer.h"
#include "linearlayer.h"
#include "logreglayer.h"
#include "minmaxlayer.h"

NeuronModel::NeuronModel(int output_size_) :
	layers(),
	se(output_size_),
	output_size(output_size_)
{

}

std::vector<double> NeuronModel::forward(const std::vector<double>& x_) {
	int batch_size = x_.size();
	std::vector<double> x = x_;
	for (int i = 0; i < layers.size(); ++i) {
		x = layers[i]->forward(x, batch_size);
	}
	return x;
}

double NeuronModel::forward(const std::vector<double>& x_, const std::vector<double>& y) {
	int batch_size = x_.size();
	std::vector<double> x = x_;
	for (int i = 0; i < layers.size(); ++i) {
		x = layers[i]->forward(x, batch_size);
	}
	std::vector<double> error = se.forward(x, y, batch_size);
	double er = 0;
	for (int i = 0; i < error.size(); ++i) {
		er += error[i];
	}
	er /= error.size();
	std::cerr << "Mean Square Error = " << er << std::endl;
	return er;
}

void NeuronModel::backward() {
	std::vector<double> dL_wrt_output = se.backward();
	for (int i = layers.size()-1; i >= 0; --i) {
		//std::cout << "Backpropagating layer " << i << std::endl;
		dL_wrt_output = layers[i]->backward(dL_wrt_output);
	}
}

void NeuronModel::step(float rate) {
	for (int i = 0; i < layers.size(); ++i) {
		layers[i]->step(rate);
	}
}

void NeuronModel::backward1() {
	std::vector<double> dL_wrt_output = se.backward();
	dL_wrt_output = layers[layers.size()-1]->backward(dL_wrt_output);
}

void NeuronModel::step1(float rate) {
	layers[layers.size()-1]->step(rate);
}


void NeuronModel::simplify() {
	for (int i = 0; i < layers.size(); ++i) {
		layers[i]->simplify();
	}
}

void NeuronModel::add_layer(std::shared_ptr<NeuronLayer> layer) {
	layers.push_back(layer);
}

void NeuronModel::pop_layer() {
	layers.pop_back();
}

void NeuronModel::save(FILE* file, bool lossy) {
	fwrite(&output_size, sizeof(int), 1, file);
	int size = layers.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		layers[i]->save(file, lossy);
	}
}

void NeuronModel::load(FILE* file) {
	fread(&output_size, sizeof(int), 1, file);
	int size;
	fread(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		NeuronLayer::NeuronLayerType type;
		fread(&type, sizeof(NeuronLayer::NeuronLayerType), 1, file);
		switch (type) {
			case NeuronLayer::NL_LIFTING:
				layers.push_back(std::shared_ptr<NeuronLayer>(new LiftingLayer()));
				break;
			case NeuronLayer::NL_LINEAR:
				layers.push_back(std::shared_ptr<NeuronLayer>(new LinearLayer()));
				break;
			case NeuronLayer::NL_LOGREG:
				layers.push_back(std::shared_ptr<NeuronLayer>(new LogregLayer()));
				break;
			case NeuronLayer::NL_MINMAX:
				layers.push_back(std::shared_ptr<NeuronLayer>(new MinMaxLayer()));
				break;
			default:
				throw "unknown layer";
				break;
		}
		layers[i]->load(file);
	}
}

void NeuronModel::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = output_size;
	array[pos++] = layers.size();
	for (int i = 0; i < layers.size(); ++i) {
		layers[i]->exportToArray(array, maxLen, pos);
	}
}

void NeuronModel::importFromArray(double* array, int len, int& pos) {
	output_size = array[pos++];
	int size = array[pos++];
	for (int i = 0; i < size; ++i) {
		NeuronLayer::NeuronLayerType type = (NeuronLayer::NeuronLayerType)array[pos++];
		switch (type) {
			case NeuronLayer::NL_LIFTING:
				layers.push_back(std::shared_ptr<NeuronLayer>(new LiftingLayer()));
				break;
			case NeuronLayer::NL_LINEAR:
				layers.push_back(std::shared_ptr<NeuronLayer>(new LinearLayer()));
				break;
			case NeuronLayer::NL_LOGREG:
				layers.push_back(std::shared_ptr<NeuronLayer>(new LogregLayer()));
				break;
			case NeuronLayer::NL_MINMAX:
				layers.push_back(std::shared_ptr<NeuronLayer>(new MinMaxLayer()));
				break;
			default:
				throw "unknown layer";
				break;
		}
		layers[i]->importFromArray(array, len, pos);
	}
}

void NeuronModel::print() {
	for (int i = 0; i < layers.size(); ++i) {
		layers[i]->print();
	}
}

void NeuronModel::set_output_size(int new_size) {
	se.set_size(new_size);
}
