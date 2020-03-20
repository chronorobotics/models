#include <assert.h>
#include <math.h>
#include "logreglayer.h"

LogregLayer::LogregLayer() {

}

LogregLayer::LogregLayer(int channels_, int sets_) :
	channels(channels_),
	sets(sets_),
	xw(),
	yw(),
	zw(),
	last_x(),
	bs(),
	dL_wrt_xw(),
	dL_wrt_yw(),
	dL_wrt_zw()
{
	xw.resize(channels*sets);
	yw.resize(channels*sets);
	zw.resize(channels*sets);
	random_init(xw, 3);
	random_init(yw, 3);
	random_init(zw, 1);

	dL_wrt_xw.resize(channels*sets);
	dL_wrt_yw.resize(channels*sets);
	dL_wrt_zw.resize(channels*sets);
}

LogregLayer::~LogregLayer() {

}

std::vector<double> LogregLayer::forward(const std::vector<double>& x, int batch_size) {
	assert(x.size() == batch_size*channels*3);
	std::vector<double> result;
	for (int b = 0; b < batch_size; ++b) {
		for (int c = 0; c < channels; ++c) {
			for (int s = 0; s < sets; ++s) {
				double foo = 0;
				assert(3*b*channels + 3*c + 2 < x.size());
				assert(c*sets + s < xw.size());
				assert(c*sets + s < yw.size());
				assert(c*sets + s < zw.size());
				foo += x[3*b*channels + 3*c] * xw[c*sets + s];
				foo += x[3*b*channels + 3*c + 1] * yw[c*sets + s];
				foo += x[3*b*channels + 3*c + 2] * zw[c*sets + s];
				result.push_back(1 / (1 + exp(foo)));
			}
		}
	}
	last_x = x;
	bs = batch_size;
	return result;
}

std::vector<double> LogregLayer::backward(const std::vector<double>& dL_wrt_output) {
	assert(dL_wrt_output.size() == bs*sets*channels);
	for (int c = 0; c < xw.size(); ++c) {
		dL_wrt_xw[c] = 0;
		dL_wrt_yw[c] = 0;
		dL_wrt_zw[c] = 0;
	}

	std::vector<double> result;
	result.resize(bs*channels*3, 0);
	for (int b = 0; b < bs; ++b) {
		for (int c = 0; c < channels; ++c) {
			for (int s = 0; s < sets; ++s) {
				double foo = 0;
				assert(3*b*channels + 3*c + 2 < last_x.size());
				assert(c*sets + s < xw.size());
				assert(c*sets + s < yw.size());
				assert(c*sets + s < zw.size());
				assert(c*sets + s < dL_wrt_xw.size());
				assert(c*sets + s < dL_wrt_yw.size());
				assert(c*sets + s < dL_wrt_zw.size());
				foo += last_x[3*b*channels + 3*c] * xw[c*sets + s];
				foo += last_x[3*b*channels + 3*c + 1] * yw[c*sets + s];
				foo += last_x[3*b*channels + 3*c + 2] * zw[c*sets + s];
				double blb = exp(foo);
				double bar = -blb / ((1 + blb)*(1 + blb)) * dL_wrt_output[b*channels*sets + c*sets + s];
				result[b*channels + 0] += bar * xw[c*sets + s];
				result[b*channels + 1] += bar * yw[c*sets + s];
				result[b*channels + 2] += bar * zw[c*sets + s];
				dL_wrt_xw[c*sets + s] += bar * last_x[3*b*channels + 3*c];
				dL_wrt_yw[c*sets + s] += bar * last_x[3*b*channels + 3*c + 1];
				dL_wrt_zw[c*sets + s] += bar * last_x[3*b*channels + 3*c + 2];
			}
		}
	}

	for (int c = 0; c < xw.size(); ++c) {
		dL_wrt_xw[c] /= bs;
		dL_wrt_yw[c] /= bs;
		dL_wrt_zw[c] /= bs;
	}

	return result;
}

void LogregLayer::step(double rate) {
	for (int c = 0; c < xw.size(); ++c) {
		assert(dL_wrt_xw[c] == dL_wrt_xw[c] && dL_wrt_yw[c] == dL_wrt_yw[c] && dL_wrt_zw[c] == dL_wrt_zw[c]);
		xw[c] -= rate * dL_wrt_xw[c];
		yw[c] -= rate * dL_wrt_yw[c];
		zw[c] -= rate * dL_wrt_zw[c];
	}
}

void LogregLayer::save(FILE* file, bool lossy) {
	NeuronLayerType type = NL_LOGREG;
	fwrite(&type, sizeof(NeuronLayerType), 1, file);
	fwrite(&channels, sizeof(int), 1, file);
	fwrite(&sets, sizeof(int), 1, file);
	int size = xw.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		fwrite(&xw[i], sizeof(double), 1, file);
		fwrite(&yw[i], sizeof(double), 1, file);
		fwrite(&zw[i], sizeof(double), 1, file);
	}
}

void LogregLayer::load(FILE* file) {
	fread(&channels, sizeof(int), 1, file);
	fread(&sets, sizeof(int), 1, file);
	int size;
	fread(&size, sizeof(int), 1, file);
	xw.resize(size);
	yw.resize(size);
	zw.resize(size);
	for (int i = 0; i < size; ++i) {
		fread(&xw[i], sizeof(double), 1, file);
		fread(&yw[i], sizeof(double), 1, file);
		fread(&zw[i], sizeof(double), 1, file);
	}
}

void LogregLayer::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = channels;
	array[pos++] = sets;
	array[pos++] = NL_LOGREG;
	array[pos++] = xw.size();
	for (int i = 0; i < xw.size(); ++i) {
		array[pos++] = xw[i];
		array[pos++] = yw[i];
		array[pos++] = zw[i];
	}
}

void LogregLayer::importFromArray(double* array, int len, int& pos) {
	channels = array[pos++];
	sets = array[pos++];
	int size = array[pos++];
	xw.resize(size);
	yw.resize(size);
	zw.resize(size);
	for (int i = 0; i < size; ++i) {
		xw[i] = array[pos++];
		yw[i] = array[pos++];
		zw[i] = array[pos++];
	}
}

void LogregLayer::print() {
	std::cout << "Logreg: ";
	for (int i = 0; i < xw.size(); ++i) {
		std::cout << "{" << xw[i] <<", " << yw[i] <<", " << zw[i] <<"}, ";
	}
	std::cout << std::endl;
}
