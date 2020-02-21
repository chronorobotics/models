#include <math.h>
#include <iostream>
#include <sstream>
#include "em_sqdist.h"

EMSqdist::EMSqdist() :
	EMCircular(),
	clusters(),
	timestamps_weight(0)
{

}

EMSqdist::EMSqdist(int cluster_count_) :
	EMCircular(cluster_count_),
	clusters(),
	timestamps_weight(0)
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

	for (unsigned int i = 0; i < timestamps.size(); ++i) {
		double s = 0;
		for (unsigned int j = 0; j < clusters.size(); ++j) {
			double a = clusters[j].density_at(timestamps[i].phase) * clusters[j].weight;
			timestamps[i].alpha[j] = a;
			s += a;
		}

		for (unsigned int j = 0; j < clusters.size(); ++j) {
			timestamps[i].alpha[j] /= s;
		}
	}
}

double EMSqdist::maximisation() {
	double shift = 0;
	std::cerr << "Performing maximisation ..." << std::endl;
	for (unsigned int i = 0; i < clusters.size(); ++i) {
		double last_xx = clusters[i].xx;
		double last_yy = clusters[i].yy;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (unsigned int j = 0; j < timestamps.size(); ++j) {
			s += timestamps[j].alpha[i];
		}
		clusters[i].weight = s / timestamps.size();

		mle(i, s, clusters[i].xx, clusters[i].yy);
		double norm = sqrt(clusters[i].xx*clusters[i].xx + clusters[i].yy*clusters[i].yy);
		if (norm > 0.999) {
			clusters[i].xx *= 0.999/norm;
			clusters[i].yy *= 0.999/norm;
		}

		double delta_xx = last_xx - clusters[i].xx;
		double delta_yy = last_yy - clusters[i].yy;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_xx*delta_xx + delta_yy*delta_yy + delta_weight*delta_weight;
	}
	return sqrt(shift);
}

void EMSqdist::train() {
	double shift = 555;
	double dl;
	int i = 5;
	do {
		double l = get_loglikelihood();
		expectation();
		shift = maximisation();
		dl = get_loglikelihood() - l;
		std::cout << "dl = " << dl << std::endl;
		if (dl < 1) {
			--i;
		}
	} while (i);
}

void EMSqdist::add_time(uint32_t time, double value) {
	timestamps.push_back(Timestamp(time, cluster_count, value));
	timestamps_weight += value;
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

double EMSqdist::get_density_at_d(double phase) const {
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(phase);
	}
	return result;
}

void EMSqdist::U(double zr, double zi, double fr, double fi, double w, double& vr, double& vi) {
	double a = zr - fr;
	double b = zi - fi;
	double c = 1 - zr*fr - zi*fi;
	double d = zr*fi - zi*fr;
	double e = w/(c*c + d*d);
	vr += e*(a*c + b*d);
	vi += e*(b*c - a*d);
}

void EMSqdist::mle(int c, double s, double& xhat, double& yhat) {
	xhat = 0;
	yhat = 0;
	int i = 0;
	while (true) {
		double vr = 0, vi = 0;
		double xhat_ = 0, yhat_ = 0;
		for (unsigned int j = 0; j < timestamps.size(); ++j) {
			U(cos(timestamps[j].phase), sin(timestamps[j].phase), xhat, yhat, timestamps[j].alpha[c], vr, vi);
		}
		U(vr/s, vi/s, -xhat, -yhat, 1, xhat_, yhat_);
		double dx = xhat - xhat_;
		double dy = yhat - yhat_;
		xhat = xhat_;
		yhat = yhat_;
		double d = sqrt(dx*dx + dy*dy);
		if (d < 1e-7) {
			return;
		}
		++i;
	}
}
