#include <math.h>
#include "neuronlayer.h"

NeuronLayer::NeuronLayer()
{

}

NeuronLayer::~NeuronLayer()
{

}

void NeuronLayer::random_init(std::vector<double>& vec, double scale) {
	for (int i = 0; i < vec.size(); ++i) {
		vec[i] = (double(rand())/RAND_MAX * 2 - 1) * scale;
	}
}

void NeuronLayer::simplify() {

}
