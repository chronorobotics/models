#include <assert.h>
#include <math.h>
#include "linearlayer.h"

LinearLayer::LinearLayer() {

}

LinearLayer::LinearLayer(int in, int out) :
	rows(out),
	cols(in+1),
	W(),
	last_x(),
	bs(),
	dL_wrt_W()
{
	W.resize(rows*cols);
	random_init(W, 1/double(cols));

	dL_wrt_W.resize(rows*cols);
}

LinearLayer::~LinearLayer() {

}

std::vector<double> LinearLayer::forward(const std::vector<double>& x, int batch_size) {
	assert(x.size() == batch_size*(cols-1));
	std::vector<double> result;
	result.resize(batch_size*rows, 0);
	int in = cols-1;
	for (int b = 0; b < batch_size; ++b) {
		for (int r = 0; r < rows; ++r) {
			for (int c = 0; c < in; ++c) {
				assert(b*rows + r < result.size());
				assert(b*in + c < x.size());
				assert(r*cols + c < W.size());
				result[b*rows + r] += x[b*in + c] * W[r*cols + c];
			}
			assert(b*rows + r < result.size());
			assert(r*cols + in < W.size());
			result[b*rows + r] += W[r*cols + in];
		}
	}

	last_x = x;
	bs = batch_size;
	return result;
}

std::vector<double> LinearLayer::backward(const std::vector<double>& dL_wrt_output) {
	assert(dL_wrt_output.size() == bs*rows);
	for (int i = 0; i < dL_wrt_W.size(); ++i) {
		dL_wrt_W[i] = 0;
	}

	std::vector<double> result;
	int in = cols-1;
	result.resize(bs*in, 0);
	for (int b = 0; b < bs; ++b) {
		for (int r = 0; r < rows; ++r) {
			for (int c = 0; c < in; ++c) {
				assert(b*in + c < result.size());
				assert(b*rows + r < dL_wrt_output.size());
				assert(r*cols + c < W.size());
				assert(r*cols + c < dL_wrt_W.size());
				assert(b*rows + r < dL_wrt_output.size());
				assert(b*in + c < last_x.size());
				result[b*in + c] += dL_wrt_output[b*rows + r] * W[r*cols + c];
				dL_wrt_W[r*cols + c] += dL_wrt_output[b*rows + r] * last_x[b*in + c];
			}
			assert(r*cols + in < dL_wrt_W.size());
			assert(b*rows + r < dL_wrt_output.size());
			dL_wrt_W[r*cols + in] += dL_wrt_output[b*rows + r];
		}
	}

	for (int i = 0; i < dL_wrt_W.size(); ++i) {
		dL_wrt_W[i] /= bs;
	}
	return result;
}

void LinearLayer::step(double rate) {
	for (int i = 0; i < W.size(); ++i) {
		assert(dL_wrt_W[i] == dL_wrt_W[i]);
		W[i] -= dL_wrt_W[i] * rate;
	}
}

void LinearLayer::save(FILE* file, bool lossy) {
	NeuronLayerType type = NL_LINEAR;
	fwrite(&type, sizeof(NeuronLayerType), 1, file);
	fwrite(&rows, sizeof(int), 1, file);
	fwrite(&cols, sizeof(int), 1, file);
	int size = W.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		fwrite(&W[i], sizeof(double), 1, file);
	}
}

void LinearLayer::load(FILE* file) {
	int size;
	fread(&rows, sizeof(int), 1, file);
	fread(&cols, sizeof(int), 1, file);
	fread(&size, sizeof(int), 1, file);
	W.resize(size);
	for (int i = 0; i < size; ++i) {
		fread(&W[i], sizeof(double), 1, file);
	}
}

void LinearLayer::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = NL_LINEAR;
	array[pos++] = rows;
	array[pos++] = cols;
	array[pos++] = W.size();
	for (int i = 0; i < W.size(); ++i) {
		array[pos++] = W[i];
	}
}

void LinearLayer::importFromArray(double* array, int len, int& pos) {
	rows = array[pos++];
	cols = array[pos++];
	int size = array[pos++];
	W.resize(size);
	for (int i = 0; i < size; ++i) {
		W[i] = array[pos++];
	}
}

void LinearLayer::print() {
	std::cout << "Linear: ";
	for (int i = 0; i < W.size(); ++i) {
		std::cout << W[i] <<", ";
	}
	std::cout << std::endl;
}

void LinearLayer::simplify() {
	for (int r = 0; r < rows; ++r) {
		double mean = 0;
		for (int c = 0; c < cols; ++c) {
			mean += fabs(W[r*cols + c]);
		}
		mean /= cols;
		double mean2 = 0;
		for (int c = 0; c < cols; ++c) {
			if (fabs(W[r*cols + c]) < mean) {
				W[r*cols + c] = 0;
			}
			mean2 += fabs(W[r*cols + c]);
		}
		mean2 /= cols;
		for (int c = 0; c < cols; ++c) {
			W[r*cols + c] *= mean/mean2;
		}
	}
}
