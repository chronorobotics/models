#include <math.h>
#include <iostream>
#include <sstream>
#include "em_unwrapped.h"

EMUnwrapped::EMUnwrapped() :
	EMCircular(),
	clusters()
{

}

EMUnwrapped::EMUnwrapped(int cluster_count_) :
	EMCircular(cluster_count_),
	clusters()
{
	double s = 0;
	cluster_count = cluster_count_;
	for (int i = 0; i < cluster_count; ++i) {
		clusters.push_back(Cluster());
		s += clusters[i].weight;
	}

	for (int i = 0; i < cluster_count; ++i) {
		clusters[i].weight /= s;
	}
}

void EMUnwrapped::expectation() {
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

double EMUnwrapped::maximisation() {
	double shift = 0;
	std::cerr << "Performing maximisation ..." << std::endl;

	for (int i = 0; i < clusters.size(); ++i) {
		double last_mu = clusters[i].mu;
		double last_sigma = clusters[i].sigma;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			s += timestamps[j].alpha[i];
		}
		clusters[i].weight = s / timestamps.size();

		double mean = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			mean += timestamps[j].phase * timestamps[j].alpha[i];
		}
		mean /= s;

		double deviation = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			double foo = timestamps[j].phase - mean;
			deviation += foo*foo * timestamps[j].alpha[i];
		}
		deviation /= s;

		clusters[i].mu = mean;
		clusters[i].sigma = sqrt(deviation);
		if (clusters[i].sigma < 0.1) clusters[i].sigma = 0.1;

		double delta_mu = last_mu - clusters[i].mu;
		double delta_sigma = last_sigma - clusters[i].sigma;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_mu*delta_mu + delta_sigma*delta_sigma + delta_weight*delta_weight;
	}
	return sqrt(shift);
}

void EMUnwrapped::train() {
	double shift = 555;
	double dl;
	int i = 5;
	do {
		std::cout << "..." << std::endl;
		for (int i = 0; i < clusters.size(); ++i) {
			if (clusters[i].sigma < 0.1) clusters[i].sigma = 0.1;
		}
		double l = get_loglikelihood();
		expectation();
		shift = maximisation();
		dl = get_loglikelihood() - l;
		std::cout << "dl = " << dl << std::endl;
		if (dl < 1) {
			--i;
		}
	} while (!isnan(shift) && i);
}

void EMUnwrapped::add_time(uint32_t time, double value) {
	timestamps.push_back(Timestamp(time, cluster_count, value));
}

EMUnwrapped::Cluster::Cluster() :
	mu(2*M_PI*float(rand()) / RAND_MAX - M_PI),
	sigma(float(rand()) / RAND_MAX/10),
	weight(float(rand()) / RAND_MAX)
{

}

EMUnwrapped::Cluster::Cluster(double mu_, double sigma_, double weight_) :
	mu(mu_),
	sigma(sigma_),
	weight()
{

}

void EMUnwrapped::Cluster::print() {
	std::cout << "(" << mu << ", " << sigma << ", " << weight << ")";
}

void EMUnwrapped::Cluster::save(FILE* file, bool lossy) {
	fwrite(&sigma, sizeof(double), 1, file);
	fwrite(&mu, sizeof(double), 1, file);
	fwrite(&weight, sizeof(double), 1, file);
}

void EMUnwrapped::Cluster::load(FILE* file) {
	fread(&sigma, sizeof(double), 1, file);
	fread(&mu, sizeof(double), 1, file);
	fread(&weight, sizeof(double), 1, file);
}

void EMUnwrapped::Cluster::importFromArray(double* array, int len, int& pos) {
	sigma = array[pos++];
	mu = array[pos++];
	weight = array[pos++];
}

void EMUnwrapped::Cluster::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = sigma;
	array[pos++] = mu;
	array[pos++] = weight;
}

double EMUnwrapped::Cluster::density_at(double value) const {
	double t = value - mu;
	return exp(-t*t / (2*sigma*sigma)) / (sigma * sqrt(2*M_PI));
}

void EMUnwrapped::print() {
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

double EMUnwrapped::get_density_at_d(double phase) const {
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(phase);
	}
	return result;
}

void EMUnwrapped::save(FILE* file, bool lossy)
{
	int size = clusters.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		clusters[i].save(file, lossy);
	}
}

void EMUnwrapped::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].load(file);
	}
}

void EMUnwrapped::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = clusters.size();
	for (int i = 0; i < clusters.size(); ++i) {
		clusters[i].exportToArray(array, maxLen, pos);
	}
}

void EMUnwrapped::importFromArray(double* array, int len, int& pos)
{
	int size = array[pos++];
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].importFromArray(array, len, pos);
	}
}
