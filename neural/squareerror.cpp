#include <assert.h>
#include <iostream>
#include "squareerror.h"

SquareError::SquareError(int size_) :
	size(size_),
	last_x(),
	last_y(),
	bs()
{

}

std::vector<double> SquareError::forward(const std::vector<double>& x, const std::vector<double>& y, int batch_size) {
	assert(x.size() == size*batch_size);
	assert(y.size() == size*batch_size);
	std::vector<double> result;
	result.reserve(batch_size);
	for (int b = 0; b < batch_size; ++b) {
		double foo = 0;
		for (int i = 0; i < size; ++i) {
			double bar = x[b*size + i] - y[b*size + i];
			foo += bar*bar;
		}
		result.push_back(foo);
	}
	bs = batch_size;
	last_x = x;
	last_y = y;
	return result;
}

std::vector<double> SquareError::backward() {
	std::vector<double> result;
	result.reserve(last_x.size());
	for (int i = 0; i < last_x.size(); ++i) {
		result.push_back(2*(last_x[i] - last_y[i]));
	}
	return result;
}
