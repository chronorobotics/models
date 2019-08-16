#include <math.h>
#include <iostream>
#include <sstream>
#include "em_multi_alternative.h"

EMMultiAlternative::EMMultiAlternative(int cluster_count_, int dimension_) :
	cluster_count(cluster_count_),
	dimension(dimension_),
	clusters(),
	measurements()
{
	double s = 0;
	for (int i = 0; i < cluster_count; ++i) {
		clusters.push_back(Cluster(dimension));
		s += clusters[i].weight;
	}

	for (int i = 0; i < cluster_count; ++i) {
		clusters[i].weight /= s;
	}
}

void EMMultiAlternative::expectation() {
	std::cerr << "Performing expectation ..." << std::endl;
	std::cout << clusters[0].mean[0] << std::endl;
	for (int i = 0; i < measurements.size(); ++i) {
		double s = 0;
		for (int j = 0; j < clusters.size(); ++j) {
			double a = clusters[j].probability_at(measurements[i].value);
			s += a;
		}

		for (int j = 0; j < clusters.size(); ++j) {
			measurements[i].alpha[j] /= s;
		}
	}
}

double EMMultiAlternative::maximisation() { //TODO
	double shift;
	std::cerr << "Performing maximisation ..." << std::endl;

	for (int i = 0; i < clusters.size(); ++i) {
		std::vector<double> last_mean = clusters[i].mean;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < measurements.size(); ++j) {
			s += measurements[j].alpha[i];
		}
		clusters[i].weight = s / measurements.size();

		std::vector<double> mean;
		mean.resize(dimension, 0);
		for (int j = 0; j < measurements.size(); ++j) {
			for (int k = 0; k < dimension; ++k) {
				if (measurements[j].value[k]) {
					mean[k] += measurements[j].alpha[i];
				}
			}
		}
		for (int k = 0; k < dimension; ++k) {
			clusters[i].mean[k] = mean[k] / s;
			double delta = last_mean[k] - clusters[i].mean[k];
			shift += delta*delta;
		}

		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_weight*delta_weight;
	}
	return sqrt(shift);
}

void EMMultiAlternative::train() {
	double shift = 555;
	std::cout << "Clustering values ..." << std::endl;
	do {
		expectation();
		shift = maximisation();
		std::cout << "shift = " << shift << std::endl;
	} while (shift > 1E-3);
}

void EMMultiAlternative::add_value(std::vector<bool> value) {
	measurements.push_back(Measurement(value, cluster_count));
}

EMMultiAlternative::Measurement::Measurement(std::vector<bool> value_, int cluster_count) :
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

EMMultiAlternative::Cluster::Cluster(int dimension_) :
	mean(),
	weight(float(rand()) / RAND_MAX)
{
	for (int i = dimension_; i; --i) {
		mean.push_back(float(rand()) / RAND_MAX);
	}
}

EMMultiAlternative::Cluster::Cluster(std::vector<double> mean_, double weight_) :
	mean(mean_),
	weight(weight_)
{

}

void EMMultiAlternative::Cluster::print() {
	std::cout << "(";
	for (int i = 0; i < mean.size(); ++i) {
		if (i) {
			std::cout << ", ";
		}
		std::cout << mean[i];
	}
	std::cout << ", w=" << weight << ")";
}

double EMMultiAlternative::Cluster::probability_at(std::vector<bool> value) const {
	double result = 1;
	for (int i = 0; i < mean.size(); ++i) {
		if (value[i]) {
			result *= mean[i];
		} else {
			result *= 1 - mean[i];
		}
	}
	return result;
}

void EMMultiAlternative::print() {
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

std::vector<double> EMMultiAlternative::get_alpha_at(std::vector<bool> value) const {
	std::vector<double> result;
	result.resize(cluster_count, 0);
	double s = 0;

	for (int i = 0; i < clusters.size(); ++i) {
		double foo = clusters[i].weight * clusters[i].probability_at(value);
		result[i] = foo;
		s += foo;
	}

	for (int i = 0; i < clusters.size(); ++i) {
		result[i] /= s;
	}

	return result;
}

std::vector<std::vector<double> > EMMultiAlternative::get_means() const {
	std::vector<std::vector<double> > result;

	for (int i = 0; i < clusters.size(); ++i) {
		result.push_back(clusters[i].mean);
	}

	return result;
}

int EMMultiAlternative::get_cluster_count() const {
	return cluster_count;
}
