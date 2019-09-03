#include <math.h>
#include <iostream>
#include <sstream>

#include "em_gauss_sqdist.h"

EMGaussSqdist::EMGaussSqdist() :
	cluster_count(),
	total_weight(0),
	clusters(),
	points()
{

}

EMGaussSqdist::EMGaussSqdist(int cluster_count_) :
	cluster_count(cluster_count_),
	total_weight(0),
	clusters(),
	points()
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

double EMGaussSqdist::time_to_phase(uint32_t time) {
	float phase = fmodf(time, 86400.0f) / 86400;
	if (phase > 0.5) {
		phase -= 1;
	}
	return phase * M_PI * 2;
}

void EMGaussSqdist::expectation() {
	std::cerr << "Performing expectation ..." << std::endl;

	for (int i = 0; i < points.size(); ++i) {
		double s = 0;
		for (int j = 0; j < clusters.size(); ++j) {
			double a = clusters[j].density_at(points[i].x, points[i].y, points[i].phase);
			points[i].alpha[j] = a;
			s += a;
		}

		for (int j = 0; j < clusters.size(); ++j) {
			points[i].alpha[j] /= s;
		}
	}
}

double EMGaussSqdist::maximisation() {
	double shift;
	std::cerr << "Performing maximisation ..." << std::endl;
	for (int i = 0; i < clusters.size(); ++i) {
		double last_ex = clusters[i].ex;
		double last_ey = clusters[i].ey;
		double last_cxx = clusters[i].cxx;
		double last_cyy = clusters[i].cyy;
		double last_cxy = clusters[i].cxy;
		double last_txx = clusters[i].txx;
		double last_tyy = clusters[i].tyy;
		double last_weight = clusters[i].weight;

		double s = 0;
		for (int j = 0; j < points.size(); ++j) {
			s += points[j].alpha[i] * points[j].weight;
		}
		clusters[i].weight = s / total_weight;

		double mean_re = 0;
		double mean_im = 0;
		double mean_x = 0;
		double mean_y = 0;
		for (int j = 0; j < points.size(); ++j) {
			mean_re += cos(points[j].phase) * points[j].alpha[i] * points[j].weight;
			mean_im += sin(points[j].phase) * points[j].alpha[i] * points[j].weight;
			mean_x += points[j].x * points[j].alpha[i] * points[j].weight;
			mean_y += points[j].y * points[j].alpha[i] * points[j].weight;
		}
		clusters[i].txx = mean_re / s;
		clusters[i].tyy = mean_im / s;
		mean_x /= s;
		mean_y /= s;
		clusters[i].ex = mean_x;
		clusters[i].ey = mean_y;
		double norm = sqrt(clusters[i].txx*clusters[i].txx + clusters[i].tyy*clusters[i].tyy);
		if (norm > 0.9) {
			clusters[i].txx *= 0.9/norm;
			clusters[i].tyy *= 0.9/norm;
		}

		double deviation_x = 0;
		double deviation_y = 0;
		double covariance = 0;
		for (int j = 0; j < points.size(); ++j) {
			double dx = points[j].x - mean_x;
			double dy = points[j].y - mean_y;
			deviation_x += dx*dx * points[j].alpha[i] * points[j].weight;
			deviation_y += dy*dy * points[j].alpha[i] * points[j].weight;
			covariance += dx*dy * points[j].alpha[i] * points[j].weight;
		}
		deviation_x /= s;
		deviation_y /= s;
		covariance /= s;
		double foo = sqrt(deviation_x*deviation_y);
		if (foo <= 0.7) {
			foo = 0.7 / foo;
			deviation_x *= foo;
			deviation_y *= foo;
			covariance *= foo;
		}

		clusters[i].det = deviation_x*deviation_y - covariance*covariance;
		clusters[i].cxx = deviation_y / clusters[i].det;
		clusters[i].cyy = deviation_x / clusters[i].det;
		clusters[i].cxy = -2*covariance / clusters[i].det;
		clusters[i].det = sqrt(4*M_PI*M_PI*clusters[i].det);

		double delta_ex = last_ex - clusters[i].ex;
		double delta_ey = last_ey - clusters[i].ey;
		double delta_cxx = last_cxx - clusters[i].cxx;
		double delta_cyy = last_cyy - clusters[i].cyy;
		double delta_cxy = last_cxy - clusters[i].cxy;
		double delta_txx = last_txx - clusters[i].txx;
		double delta_tyy = last_tyy - clusters[i].tyy;
		double delta_weight = last_weight - clusters[i].weight;
		shift += delta_ex*delta_ex + delta_ey*delta_ey + delta_cxx*delta_cxx + delta_cyy*delta_cyy +
						 delta_cxy*delta_cxy + delta_txx*delta_txx + delta_tyy*delta_tyy + delta_weight*delta_weight;
	}
	return sqrt(shift);
}

void EMGaussSqdist::train() {
	double shift = 555;
	do {
		expectation();
		shift = maximisation();
		std::cerr << "shift = " << shift << std::endl;
	} while (shift > 1E-3);
}

void EMGaussSqdist::add_value(double x, double y, uint32_t time, double weight) {
	points.push_back(Point(x, y, time, cluster_count, weight));
	total_weight += weight;
}

EMGaussSqdist::Cluster::Cluster() :
	ex(2*float(rand()) / RAND_MAX - 1),
	ey(2*float(rand()) / RAND_MAX - 1),
	cxx(),
	cyy(),
	cxy(),
	det(),
	txx(),
	tyy(),
	weight(float(rand()) / RAND_MAX)
{
	double r = double(rand())/RAND_MAX;
	double f = 2*M_PI*double(rand())/RAND_MAX;
	txx = r*cos(f);
	tyy = r*sin(f);
	double a = float(rand()) / RAND_MAX + 1;
	double b = float(rand()) / RAND_MAX + 1;
	double c = float(rand()) / RAND_MAX * sqrt(a*b);
	det = a*b - c*c;
	cxx = b / det;
	cyy = a / det;
	cxy = -2*c / det;
	det = sqrt(4*M_PI*M_PI*det);
}

EMGaussSqdist::Cluster::Cluster(double ex_, double ey_, double cxx_, double cyy_,
																double cxy_, double det_, double txx_, double tyy_,
																double weight_) :
	ex(ex_),
	ey(ey_),
	cxx(cxx_),
	cyy(cyy_),
	cxy(cxy_),
	det(det_),
	txx(txx_),
	tyy(tyy_),
	weight(weight_)
{

}

void EMGaussSqdist::Cluster::print() {
	std::cout << "(" << ex << ", " << ey << ", " << cxx << ", " << cyy << ", " << cxy << ", " << txx << ", " << tyy << weight << ")";
}

void EMGaussSqdist::Cluster::save(FILE* file, bool lossy) {
	fwrite(&ex, sizeof(double), 1, file);
	fwrite(&ey, sizeof(double), 1, file);
	fwrite(&cxx, sizeof(double), 1, file);
	fwrite(&cxy, sizeof(double), 1, file);
	fwrite(&cyy, sizeof(double), 1, file);
	fwrite(&det, sizeof(double), 1, file);
	fwrite(&txx, sizeof(double), 1, file);
	fwrite(&tyy, sizeof(double), 1, file);
	fwrite(&weight, sizeof(double), 1, file);
}

void EMGaussSqdist::Cluster::load(FILE* file) {
	fread(&ex, sizeof(double), 1, file);
	fread(&ey, sizeof(double), 1, file);
	fread(&cxx, sizeof(double), 1, file);
	fread(&cxy, sizeof(double), 1, file);
	fread(&cyy, sizeof(double), 1, file);
	fread(&det, sizeof(double), 1, file);
	fread(&txx, sizeof(double), 1, file);
	fread(&tyy, sizeof(double), 1, file);
	fread(&weight, sizeof(double), 1, file);
}

void EMGaussSqdist::Cluster::importFromArray(double* array, int len, int& pos) {
	ex = array[pos++];
	ey = array[pos++];
	cxx = array[pos++];
	cyy = array[pos++];
	cxy = array[pos++];
	det = array[pos++];
	txx = array[pos++];
	tyy = array[pos++];
	weight = array[pos++];
}

void EMGaussSqdist::Cluster::exportToArray(double* array, int maxLen, int& pos) {
	array[pos++] = ex;
	array[pos++] = ey;
	array[pos++] = cxx;
	array[pos++] = cyy;
	array[pos++] = cxy;
	array[pos++] = det;
	array[pos++] = txx;
	array[pos++] = tyy;
	array[pos++] = weight;
}

double EMGaussSqdist::Cluster::density_at(double x, double y, double phase) const {
	double dtx = txx - cos(phase);
	double dty = tyy - sin(phase);
	double sqdist = (1 - txx*txx - tyy*tyy) / ((dtx*dtx + dty*dty) * 2*M_PI);
	double dx = ex - x;
	double dy = ey - y;
	double gauss = exp(-(cxx*dx*dx + cyy*dy*dy + cxy*dx*dy)/2) / det;
	return sqdist * gauss;
}

void EMGaussSqdist::save(FILE* file, bool lossy)
{
	int size = clusters.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		clusters[i].save(file, lossy);
	}
}

void EMGaussSqdist::load(FILE* file)
{
	int size;
	fread(&size, sizeof(int), 1, file);
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].load(file);
	}
}

void EMGaussSqdist::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = clusters.size();
	for (int i = 0; i < clusters.size(); ++i) {
		clusters[i].exportToArray(array, maxLen, pos);
	}
}

void EMGaussSqdist::importFromArray(double* array, int len, int& pos)
{
	int size = array[pos++];
	clusters.resize(size);
	for (int i = 0; i < size; ++i) {
		clusters[i].importFromArray(array, len, pos);
	}
}

void EMGaussSqdist::print() {
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

double EMGaussSqdist::get_density_at(double x, double y, uint32_t time) const {
	double phase = time_to_phase(time);
	double result = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		result += clusters[i].weight * clusters[i].density_at(x, y, phase);
	}
	return result;
}

EMGaussSqdist::Point::Point(double x_, double y_, uint32_t time, int cluster_count, double weight_) :
	x(x_),
	y(y_),
	phase(time_to_phase(time)),
	alpha(),
	weight(weight_)
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
