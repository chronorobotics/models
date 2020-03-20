#include <algorithm>
#include <assert.h>
#include "minmaxlayer.h"

MinMaxLayer::MinMaxLayer() {

}

MinMaxLayer::MinMaxLayer(int size_) :
	size(size_),
	last_x(),
	bs()
{

}

MinMaxLayer::~MinMaxLayer() {

}

std::vector<double> MinMaxLayer::forward(const std::vector<double>& x, int batch_size) {
	assert(x.size() == size*batch_size);
	std::vector<double> result;
	for (int b = 0; b < batch_size; ++b) {
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				assert(b*size + i < x.size());
				assert(b*size + j < x.size());
				if (i < j) {
					result.push_back(std::max(x[b*size + i], x[b*size + j]));
				} else if (i == j) {
					result.push_back(x[b*size + i]);
				} else {
					result.push_back(std::min(x[b*size + i], x[b*size + j]));
				}
			}
		}
	}
	last_x = x;
	bs = batch_size;
	return result;
}

std::vector<double> MinMaxLayer::backward(const std::vector<double>& dL_wrt_output) {
	assert(dL_wrt_output.size() == size*size*bs);
	int s2 = size*size;
	std::vector<double> result;
	result.resize(bs*size, 0);
	for (int b = 0; b < bs; ++b) {
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				assert(b*s2 + i*size + j < dL_wrt_output.size());
				assert(b*size + i < last_x.size());
				assert(b*size + j < last_x.size());
				assert(b*size + i < result.size());
				assert(b*size + j < result.size());
				double foo = dL_wrt_output[b*s2 + i*size + j];
				double xi = last_x[b*size + i];
				double xj = last_x[b*size + j];
				if (i < j) {
					if (xi < xj) {
						result[b*size + i] += foo;
					} else {
						result[b*size + j] += foo;
					}
				} else if (i == j) {
					result[b*size + i] += foo;
				} else {
					if (xi > xj) {
						result[b*size + i] += foo;
					} else {
						result[b*size + j] += foo;
					}
				}
			}
		}
	}
	return result;
}

void MinMaxLayer::step(double rate) {

}

void MinMaxLayer::save(FILE* file, bool lossy) {
	NeuronLayerType type = NL_MINMAX;
	fwrite(&type, sizeof(NeuronLayerType), 1, file);
	fwrite(&size, sizeof(int), 1, file);
}

void MinMaxLayer::load(FILE* file) {
	fread(&size, sizeof(int), 1, file);
}

void MinMaxLayer::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = NL_MINMAX;
	array[pos++] = size;
}

void MinMaxLayer::importFromArray(double* array, int len, int& pos) {
	size = array[pos++];
}

void MinMaxLayer::print() {
	std::cout << "Minmax" << std::endl;
}
