#include <assert.h>
#include <math.h>

#include "liftinglayer.h"

LiftingLayer::LiftingLayer() {

}

LiftingLayer::LiftingLayer(std::vector<double> periods_, double t_min, double t_max) :
	periods(periods_),
	scale(2/(t_max-t_min)),
	offset((t_max+t_min)/2),
	last_x(),
	bs()
{

}

LiftingLayer::~LiftingLayer() {

}

std::vector<double> LiftingLayer::forward(const std::vector<double>& x, int batch_size) {
	assert(x.size() == batch_size);
	std::vector<double> result;
	result.reserve(3*batch_size*(periods.size()+1));
	for (int b = 0; b < batch_size; ++b) {
		for (int p = 0; p < periods.size(); ++p) {
			double phase = fmodf(x[b], periods[p]) / periods[p] * 2 * M_PI;
			result.push_back(cos(phase));
			result.push_back(sin(phase));
			result.push_back(1);
		}
		double x_ = (x[b] - offset)*scale;
		result.push_back(x_);
		result.push_back(x_*x_);
		result.push_back(1);
	}
	last_x = x;
	bs = batch_size;
	return result;
}

std::vector<double> LiftingLayer::backward(const std::vector<double>& dL_wrt_output) {
	assert(dL_wrt_output.size() == 3*bs*(1+periods.size()));
	std::vector<double> result;
	result.resize(bs, 0);
	int ps1 = periods.size()+1;
	int ps = periods.size();
	for (int b = 0; b < bs; ++b) {
		double foo = 0;
		for (int p = 0; p < ps; ++p) {
			double phase = fmodf(last_x[b], periods[p]) / periods[p] * 2 * M_PI;
			foo += dL_wrt_output[3*b*ps1 + 3*p] * (-2*M_PI/periods[p]*sin(phase));
			foo += dL_wrt_output[3*b*ps1 + 3*p + 1] * (2*M_PI/periods[p]*cos(phase));
		}
		foo += dL_wrt_output[3*b*ps1 + 3*ps] * scale;
		foo += dL_wrt_output[3*b*ps1 + 3*ps + 1] * 2*scale*scale*(last_x[b] - offset);
		result[b] = foo;
	}
	return result;
}

void LiftingLayer::step(double rate) {

}

void LiftingLayer::save(FILE* file, bool lossy) {
	NeuronLayerType type = NL_LIFTING;
	fwrite(&type, sizeof(NeuronLayerType), 1, file);
	fwrite(&scale, sizeof(double), 1, file);
	fwrite(&offset, sizeof(double), 1, file);
	int size = periods.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		fwrite(&periods[i], sizeof(double), 1, file);
	}
}

void LiftingLayer::load(FILE* file) {
	int size;
	fread(&scale, sizeof(double), 1, file);
	fread(&offset, sizeof(double), 1, file);
	fread(&size, sizeof(int), 1, file);
	periods.resize(size);
	for (int i = 0; i < size; ++i) {
		fread(&periods[i], sizeof(double), 1, file);
	}
}

void LiftingLayer::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = NL_LIFTING;
	array[pos++] = scale;
	array[pos++] = offset;
	array[pos++] = periods.size();
	for (int i = 0; i < periods.size(); ++i) {
		array[pos++] = periods[i];
	}
}

void LiftingLayer::importFromArray(double* array, int len, int& pos) {
	scale = array[pos++];
	offset = array[pos++];
	int size = array[pos++];
	periods.resize(size);
	for (int i = 0; i < size; ++i) {
		periods[i] = array[pos++];
	}
}

void LiftingLayer::print() {
	std::cout << "Lifting: ";
	for (int i = 0; i < periods.size(); ++i) {
		std::cout << periods[i] <<", ";
	}
	std::cout << std::endl;
}
