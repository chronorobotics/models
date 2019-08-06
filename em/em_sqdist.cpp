#include <math.h>
#include <iostream>
#include <sstream>
#include "em_sqdist.h"

EMSqdist::EMSqdist(int cluster_count_) :
	cluster_count(cluster_count_),
	clusters(),
	timestamps()
{
	double s = 0;
	for (int i = 0; i < cluster_count; ++i) {
		clusters.push_back(Cluster());
		s += clusters[i].weight;
		double r = double(rand())/RAND_MAX;
		double f = 2*M_PI*double(rand())/RAND_MAX;
		clusters[i].xx = r*cos(f);
		clusters[i].yy = r*sin(f);
	}

	for (int i = 0; i < cluster_count; ++i) {
		clusters[i].weight /= s;
	}
}

void EMSqdist::expectation() {
	std::cerr << "Performing expectation ..." << std::endl;

	for (int i = 0; i < timestamps.size(); ++i) {
		double s = 0;
		for (int j = 0; j < clusters.size(); ++j) {
			double a = clusters[j].density_at(timestamps[i].phase);
			timestamps[i].alpha[j] = a;
			s += a;
		}

		for (int j = 0; j < clusters.size(); ++j) {
			timestamps[i].alpha[j] /= s;
		}
	}
}

double EMSqdist::maximisation() {
	double shift;
	std::cerr << "Performing maximisation ..." << std::endl;

	for (int i = 0; i < clusters.size(); ++i) {
		double last_xx = clusters[i].xx;
		double last_yy = clusters[i].yy;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			s += timestamps[j].alpha[i];
		}
		clusters[i].weight = s / timestamps.size();

		double mean_re = 0;
		double mean_im = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			mean_re += cos(timestamps[j].phase) * timestamps[j].alpha[i];
			mean_im += sin(timestamps[j].phase) * timestamps[j].alpha[i];
		}
		clusters[i].xx = mean_re / s;
		clusters[i].yy = mean_im / s;
		double norm = sqrt(clusters[i].xx*clusters[i].xx + clusters[i].yy*clusters[i].yy);
		if (norm > 0.99) {
			clusters[i].xx *= 0.99/norm;
			clusters[i].yy *= 0.99/norm;
		}

		double delta_xx = last_xx - clusters[i].xx;
		double delta_yy = last_yy - clusters[i].yy;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_xx*delta_xx + delta_yy*delta_yy + delta_weight*delta_weight;
	}

	return sqrt(shift);
}

double EMSqdist::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 86400.0f) / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

void EMSqdist::train() {
	double shift = 555;
	do {
		expectation();
		shift = maximisation();
		std::cout << "shift = " << shift << std::endl;
	} while (isnan(shift) || shift > 1E-3);
}

void EMSqdist::add_time(uint32_t time) {
	timestamps.push_back(Timestamp(time, cluster_count));
}

EMSqdist::Timestamp::Timestamp(uint32_t time, int cluster_count) :
	phase(time_to_phase(time)),
	alpha()
{
	double s = 0;
	double r;

	for (int i = cluster_count; i; --i) {
		r = float(rand()) / RAND_MAX;
		s += r;
		alpha.push_back(r);
	}

	for (int i = cluster_count - 1; i >= 0; --i) {
		alpha[i] /= s;
	}
}

EMSqdist::Cluster::Cluster() :
	xx(),
	yy(),
	weight(float(rand()) / RAND_MAX)
{
	double r = double(rand())/RAND_MAX;
	double f = 2*M_PI*double(rand())/RAND_MAX;
	xx = r*cos(f);
	yy = r*sin(f);
}

EMSqdist::Cluster::Cluster(double xx_, double yy_, double weight_) :
	xx(xx_),
	yy(yy_),
	weight()
{

}

void EMSqdist::Cluster::print() {
	std::cout << "(" << xx << ", " << yy << ", " << weight << ")";
}

void EMSqdist::Cluster::save(FILE* file, bool lossy) {
	fwrite(&xx, sizeof(double), 1, file);
	fwrite(&yy, sizeof(double), 1, file);
	fwrite(&weight, sizeof(double), 1, file);
}

void EMSqdist::Cluster::load(FILE* file) {
	fread(&xx, sizeof(double), 1, file);
	fread(&yy, sizeof(double), 1, file);
	fread(&weight, sizeof(double), 1, file);
}

void EMSqdist::Cluster::importFromArray(double* array, int len, int& pos) {
	xx = array[pos++];
	yy = array[pos++];
	weight = array[pos++];
}

void EMSqdist::Cluster::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = xx;
	array[pos++] = yy;
	array[pos++] = weight;
}

double EMSqdist::Cluster::density_at(double phase) const {
	double dx = xx - cos(phase);
	double dy = yy - sin(phase);
	return (1 - xx*xx - yy*yy) / ((dx*dx + dy*dy) * 2*M_PI);
}

void EMSqdist::save(FILE* file, bool lossy)
{
	int size = clusters.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		clusters[i].save(file, lossy);
	}
}

void EMSqdist::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].load(file);
	}
}

void EMSqdist::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = clusters.size();
	for (int i = 0; i < clusters.size(); ++i) {
		clusters[i].exportToArray(array, maxLen, pos);
	}
}

void EMSqdist::importFromArray(double* array, int len, int& pos)
{
	int size = array[pos++];
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].importFromArray(array, len, pos);
	}
}

void EMSqdist::print() {
	std::cout << "[";
	double sum_w = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		if (i) {
			std::cout << ", ";
		}
		clusters[i].print();
		sum_w += clusters[i].weight;
	}
	std::cout << "] (w = " << sum_w << ")" << std::endl;
}

double EMSqdist::get_density_at(uint32_t time) {
	double phase = time_to_phase(time);
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(phase);
	}
	return result;
}
