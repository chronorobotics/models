#include <math.h>
#include <iostream>
#include <sstream>
#include <ctime>

#include "em_elliptic_sqdist.h"

EMEllipticSqdist::EMEllipticSqdist() :
	EMCircular(),
	clusters(),
	timestamps_weight(0),
	test()
{

}

EMEllipticSqdist::EMEllipticSqdist(int cluster_count_, std::string filename) :
	EMCircular(cluster_count_),
	clusters(),
	timestamps_weight(0),
	test(filename)
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

EMEllipticSqdist::~EMEllipticSqdist() {
	test.close();
}

void EMEllipticSqdist::expectation() {
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

double EMEllipticSqdist::maximisation() {
	double shift = 0;
	std::cerr << "Performing maximisation ..." << std::endl;
	for (int i = 0; i < clusters.size(); ++i) {
		double last_xx = clusters[i].xx * cos(clusters[i].t0);
		double last_yy = clusters[i].xx * sin(clusters[i].t0);
		double last_aa = clusters[i].aa;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			s += timestamps[j].alpha[i] * timestamps[j].weight;
		}
		clusters[i].weight = s / timestamps_weight;

		double m1_re = 0;
		double m1_im = 0;
		double m2_re = 0;
		double m2_im = 0;
		for (int j = 0; j < timestamps.size(); ++j) {
			double foo = timestamps[j].alpha[i] * timestamps[j].weight;
			m1_re += cos(timestamps[j].phase) * foo;
			m1_im += sin(timestamps[j].phase) * foo;
			m2_re += cos(2*timestamps[j].phase) * foo;
			m2_im += sin(2*timestamps[j].phase) * foo;
		}
		m1_re /= s;
		m1_im /= s;
		m2_re /= s;
		m2_im /= s;

		double x2 = m1_re*m1_re + m1_im*m1_im;

		clusters[i].xx = sqrt(x2);
		clusters[i].t0 = atan2(m1_im, m1_re);

		double m1r2 = m1_re*m1_re - m1_im*m1_im;
		double m1i2 = 2*m1_re*m1_im;
		double m2 = sqrt(m2_re*m2_re + m2_im*m2_im);
		if (fabs(m1r2-m2_re) + fabs(m1i2-m2_im) > fabs(m1r2+m2_re) + fabs(m1i2+m2_im)) {
			m2 = -m2;
		}

		clusters[i].aa = (2*x2 - m2 - 1) / (m2 - 1);

		if (clusters[i].xx > 0.99) {
			clusters[i].xx = 0.99;
		}
		if (clusters[i].aa > 1) {
			clusters[i].aa = 1;
		}
		if (clusters[i].aa < 0.05) {
			clusters[i].aa = 0.05;
		}
		//clusters[i].aa = sqrt(1 - x2);

		test << m1_re << " " << m1_im << " " << clusters[i].aa << " " << clusters[i].weight << " ";

		double delta_xx = last_xx - clusters[i].xx * cos(clusters[i].t0);
		double delta_yy = last_yy - clusters[i].xx * sin(clusters[i].t0);
		double delta_aa = last_aa - clusters[i].aa;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_xx*delta_xx + delta_yy*delta_yy + delta_aa*delta_aa + delta_weight*delta_weight;
		//shift = std::max(shift, std::max(delta_aa, std::max(delta_weight, std::max(delta_xx, delta_yy))));
	}
	test << std::endl;
	return sqrt(shift);
}

void EMEllipticSqdist::train() {
	double shift = 555;
	int i = 0;
	do {
		expectation();
		shift = maximisation();
		std::cout << "shift = " << shift << " i = " << i << std::endl;
		++i;
	} while (shift > 1E-3 /*&& i < 1000*/);
}

void EMEllipticSqdist::add_time(uint32_t time, double value) {
	timestamps.push_back(Timestamp(time, cluster_count, value));
	timestamps_weight += value;
}

EMEllipticSqdist::Cluster::Cluster() :
	xx(double(rand())/RAND_MAX),
	t0(2*M_PI*double(rand())/RAND_MAX),
	aa(float(rand()) / RAND_MAX /*sqrt(1 - xx*xx)*/),
	weight(float(rand()) / RAND_MAX)
{

}

EMEllipticSqdist::Cluster::Cluster(double xx_, double t0_, double aa_, double weight_) :
	xx(xx_),
	t0(t0_),
	aa(aa_),
	weight()
{

}

void EMEllipticSqdist::Cluster::print() {
	std::cout << "(" << xx << ", " << t0 << ", " << aa << "," << weight << ")";
}

void EMEllipticSqdist::Cluster::save(FILE* file, bool lossy) {
	fwrite(&xx, sizeof(double), 1, file);
	fwrite(&t0, sizeof(double), 1, file);
	fwrite(&aa, sizeof(double), 1, file);
	fwrite(&weight, sizeof(double), 1, file);
}

void EMEllipticSqdist::Cluster::load(FILE* file) {
	fread(&xx, sizeof(double), 1, file);
	fread(&t0, sizeof(double), 1, file);
	fread(&aa, sizeof(double), 1, file);
	fread(&weight, sizeof(double), 1, file);
}

void EMEllipticSqdist::Cluster::importFromArray(double* array, int len, int& pos) {
	xx = array[pos++];
	t0 = array[pos++];
	aa = array[pos++];
	weight = array[pos++];
}

void EMEllipticSqdist::Cluster::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = xx;
	array[pos++] = t0;
	array[pos++] = aa;
	array[pos++] = weight;
}

double EMEllipticSqdist::Cluster::density_at(double phase) const {
	/*double dx = xx - cos(phase);
	double dy = yy - sin(phase);
	return (1 - xx*xx - yy*yy) / ((dx*dx + dy*dy) * 2*M_PI);*/
	phase -= t0;
	double dx = cos(phase) - xx;
	double dy = sin(phase) * aa;
	return aa*(1 - xx*xx) / (2 * M_PI * (dx*dx + dy*dy));
}

void EMEllipticSqdist::save(FILE* file, bool lossy)
{
	int size = clusters.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		clusters[i].save(file, lossy);
	}
}

void EMEllipticSqdist::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].load(file);
	}
}

void EMEllipticSqdist::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = clusters.size();
	for (int i = 0; i < clusters.size(); ++i) {
		clusters[i].exportToArray(array, maxLen, pos);
	}
}

void EMEllipticSqdist::importFromArray(double* array, int len, int& pos)
{
	int size = array[pos++];
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].importFromArray(array, len, pos);
	}
}

void EMEllipticSqdist::print() {
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

double EMEllipticSqdist::get_density_at(uint32_t time) {
	double phase = time_to_phase(time);
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(phase);
	}
	return result;
}
