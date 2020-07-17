#include <math.h>
#include <iostream>
#include <sstream>
#include "em_gaussian.h"

EMGaussian::EMGaussian(int cluster_count_) :
	cluster_count(cluster_count_),
	clusters(),
	measurements()
{
	double s = 0;
	for (int i = 0; i < cluster_count; ++i) {
		clusters.push_back(Cluster());
		s += clusters[i].weight;
	}

	for (int i = 0; i < cluster_count; ++i) {
		clusters[i].weight /= s;
	}
}

void EMGaussian::expectation() {
	std::cerr << "Performing expectation ..." << std::endl;

	for (int i = 0; i < measurements.size(); ++i) {
		double s = 0;
		for (int j = 0; j < clusters.size(); ++j) {
			double a = clusters[j].density_at(measurements[i].value);
			measurements[i].alpha[j] = a;
			s += a;
		}

		for (int j = 0; j < clusters.size(); ++j) {
			measurements[i].alpha[j] /= s;
		}
	}
}

double EMGaussian::maximisation() {
	double shift = 0;
	std::cerr << "Performing maximisation ..." << std::endl;

	for (int i = 0; i < clusters.size(); ++i) {
		double last_mu = clusters[i].mu;
		double last_sigma = clusters[i].sigma;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < measurements.size(); ++j) {
			s += measurements[j].alpha[i];
		}
		clusters[i].weight = s / measurements.size();

		double mean = 0;
		for (int j = 0; j < measurements.size(); ++j) {
			mean += measurements[j].value * measurements[j].alpha[i];
		}
		mean /= s;

		double deviation = 0;
		for (int j = 0; j < measurements.size(); ++j) {
			double foo = measurements[j].value - mean;
			deviation += foo*foo * measurements[j].alpha[i];
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

void EMGaussian::train() {
	double shift = 555;
	std::cout << "Clustering values ..." << std::endl;
	do {
		expectation();
		shift = maximisation();
		std::cout << "shift = " << shift << std::endl;
	} while (shift > 1E-3);
}

void EMGaussian::add_value(double value) {
	measurements.push_back(Measurement(value, cluster_count));
}

EMGaussian::Measurement::Measurement(double value_, int cluster_count) :
	value(value_),
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

EMGaussian::Cluster::Cluster() :
	mu(30*float(rand()) / RAND_MAX),
	sigma(float(rand()) / RAND_MAX),
	weight(float(rand()) / RAND_MAX)
{

}

EMGaussian::Cluster::Cluster(double mu_, double sigma_, double weight_) :
	mu(mu_),
	sigma(sigma_),
	weight()
{

}

void EMGaussian::Cluster::print() {
	std::cout << "(" << mu << ", " << sigma << ", " << weight << ")";
}

double EMGaussian::Cluster::density_at(double value) const {
	double t = value - mu;
	return exp(-t*t / (2*sigma*sigma)) / (sigma * sqrt(2*M_PI));
}

void EMGaussian::print() {
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

std::vector<double> EMGaussian::get_alpha_at(double value) const {
	std::vector<double> result;
	result.resize(cluster_count, 0);
	double s = 0;

	for (int i = 0; i < clusters.size(); ++i) {
		double foo = clusters[i].weight * clusters[i].density_at(value);
		result[i] = foo;
		s += foo;
	}

	for (int i = 0; i < clusters.size(); ++i) {
		result[i] /= s;
	}

	return result;
}

std::vector<double> EMGaussian::get_means() const {
	std::vector<double> result;

	for (int i = 0; i < clusters.size(); ++i) {
		result.push_back(clusters[i].mu);
	}

	return result;
}

int EMGaussian::get_cluster_count() const {
	return cluster_count;
}
