#include <math.h>
#include <stdio.h>
#include <iostream>

#include "logreg.hpp"

Logreg::Logreg() :
	my_weight(),
	w(0, 0, 0),
	points(),
	total_weight(),
	period(),
	scale(),
	offset(),
	parabolic(),
	unlimited()
{

}

Logreg::Logreg(double period_) :
	my_weight(1),
	w(1, 0.3, 0),
	points(),
	total_weight(0),
	period(period_),
	scale(),
	offset(),
	parabolic(false),
	unlimited(false)
{

}

Logreg::Logreg(double scale_, double offset_) :
	my_weight(1),
	w(1, 0.3, 0),
	points(),
	total_weight(0),
	period(),
	scale(scale_),
	offset(offset_),
	parabolic(true),
	unlimited(false)
{

}

double Logreg::train() {
	Vec3D w_new(1, 0, 0);
	Vec3D w_prev(1, 0.3, 0);
	double step_size = 1.0;
	double E = loss(w);
	Vec3D g = gradient(w);
	double dist = 100;
	while (dist > 0.001) {
		w_new = w_prev - g * step_size;
		double E_new = loss(w_new);
		if (E_new != E_new) {
			w_new = w_prev;
			unlimited = true;
			std::cout << "unlimited" << std::endl;
			break;
		}
		Vec3D g_new = gradient(w_new);

		if (E_new <= E) {
			dist = (w_new - w_prev).norm();
			w_prev = w_new;
			E = E_new;
			g = g_new;
			step_size *= 2;
		} else {
			step_size /= 2;
		}
	}

	w = w_new;
	//std::cout << "w = ["<< w.x << ", " << w.y << ", " << w.z << "]" << std::endl;
	/*double my_loss = 0;
	for (int i = 0; i < points.size(); ++i) {
		my_loss += ((1/(1 + exp(-(w*points[i].x))) > 0.5 ? 1 : -1) == points[i].y ? 0 : 1) * points[i].weight;
	}*/
	double my_loss = loss(w);
	points.clear();

	return my_loss;
}

void Logreg::add_measurement(double time, bool positive, float weight) {
	points.push_back(Measurement(time_to_point(time), positive ? 1 : -1, weight));
	total_weight += weight;
}

double Logreg::loss(Vec3D at) const {
	double result = 0;
	for (int i = points.size(); i >= 0; --i) {
		result += log(1 + exp(-(at*points[i].x)*points[i].y)) * points[i].weight;
	}
	return result / total_weight;
}

Logreg::Vec3D Logreg::gradient(Vec3D at) const {
	Vec3D result(0, 0, 0);
	for (int i = points.size(); i >= 0; --i) {
		result += points[i].x * (points[i].y/(1 + exp((at*points[i].x)*points[i].y))) * points[i].weight;
	}
	result.x /= -total_weight;
	result.y /= -total_weight;
	result.z /= -total_weight;
	return result;
}

double Logreg::evaluate(double time) const {
	double result = 1/(1 + exp(-(w*time_to_point(time))));
	if (unlimited) {
		result = result > 0.5 ? 1 : 0;
	}
	return result*2 - 1;
}

Logreg::Vec3D Logreg::time_to_point(double time) const {
	if (parabolic) {
		double foo = (time - offset)/scale;
		return Vec3D(foo, foo*foo, 1);
	} else {
		double phase = fmod(time, period) / period * 2*M_PI;
		return Vec3D(cos(phase), sin(phase), 1);
	}
}

void Logreg::save(FILE* file, bool lossy)
{
	fwrite(&my_weight, sizeof(double), 1, file);
	fwrite(&w, sizeof(Vec3D), 1, file);
	fwrite(&unlimited, sizeof(bool), 1, file);
	fwrite(&parabolic, sizeof(bool), 1, file);
	if (parabolic) {
		fwrite(&scale, sizeof(double), 1, file);
		fwrite(&offset, sizeof(double), 1, file);
	} else {
		fwrite(&period, sizeof(double), 1, file);
	}
}

void Logreg::load(FILE* file)
{
	fread(&my_weight, sizeof(double), 1, file);
	fread(&w, sizeof(Vec3D), 1, file);
	fread(&unlimited, sizeof(bool), 1, file);
	fread(&parabolic, sizeof(bool), 1, file);
	if (parabolic) {
		fread(&scale, sizeof(double), 1, file);
		fread(&offset, sizeof(double), 1, file);
	} else {
		fread(&period, sizeof(double), 1, file);
	}
}

void Logreg::exportToArray(double* array, int maxLen, int& pos)
{
	array[pos++] = my_weight;
	array[pos++] = w.x;
	array[pos++] = w.y;
	array[pos++] = w.z;
	array[pos++] = unlimited;
	array[pos++] = parabolic;
	if (parabolic) {
		array[pos++] = scale;
		array[pos++] = offset;
	} else {
		array[pos++] = period;
	}
}

void Logreg::importFromArray(double* array, int len, int& pos)
{
	my_weight = array[pos++];
	w.x = array[pos++];
	w.y = array[pos++];
	w.z = array[pos++];
}

void Logreg::print() {
	if (parabolic) {
		std::cout << "{parabolic, x=" << w.x << ", y=" << w.y << ", z=" << w.z << "}";
	} else {
		std::cout << "{T=" << period << ", x=" << w.x << ", y=" << w.y << ", z=" << w.z << "}";
	}
}
