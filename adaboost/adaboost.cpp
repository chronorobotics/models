#include <limits>
#include <iostream>
#include <fstream>
#include <ostream>
#include <stdio.h>
#include "adaboost.h"

Adaboost::Adaboost() :
	time_samples(),
	periods(),
	weak_classifiers(),
	min_time(std::numeric_limits<long int>::max()),
	max_time(std::numeric_limits<long int>::min()),
	sum_alpha(0)
{

}

void Adaboost::add_time(long time, bool positive) {
	if (time < min_time) {
		min_time = time;
	}
	if (time > max_time) {
		max_time = time;
	}
	time_samples.push_back(TimeSample(time, 1, positive));
}

void Adaboost::add_period(double period) {
	periods.push_back(period);
}

void Adaboost::train() {
	for (int i = 0; i < time_samples.size(); ++i) {
		time_samples[i].weight = 1/double(time_samples.size());
	}

	for (int i = 0; i < 20; ++i) {
		adaboost_iteration();
	}

	std::ofstream myfile("2.txt");
	for (int i = 0; i < time_samples.size(); ++i) {
		myfile << time_samples[i].weight << std::endl;
	}
	myfile.close();

	time_samples.clear();
}

void Adaboost::adaboost_iteration() {
	//Parabolic classifier
	Logreg wl_best((min_time + max_time)/2, 2/(max_time - min_time));
	for (int j = 0; j < time_samples.size(); ++j) {
		//std::cout << time_samples[j].time << " " << time_samples[j].positive << " " << time_samples[j].weight << std::endl;
		wl_best.add_measurement(time_samples[j].time, time_samples[j].positive, time_samples[j].weight);
	}
	double min_loss = wl_best.train();

	//Circular classifiers
	for (int i = 0; i < periods.size(); ++i) {
		Logreg wl(periods[i]);
		for (int j = 0; j < time_samples.size(); ++j) {
			wl.add_measurement(time_samples[j].time, time_samples[j].positive, time_samples[j].weight);
		}
		double loss = wl.train();
		if (loss < min_loss) {
			wl_best = wl;
			min_loss = loss;
		}
	}

	//Calculate new weights
	double mu = 0;
	for (int j = 0; j < time_samples.size(); ++j) {
		mu += time_samples[j].weight * (time_samples[j].positive ? 1 : -1) * wl_best.evaluate(time_samples[j].time);
	}
	wl_best.my_weight = log((1+mu)/(1-mu))/2;
	sum_alpha += wl_best.my_weight;

	for (int j = 0; j < time_samples.size(); ++j) {
		time_samples[j].weight *= (1 - mu*(time_samples[j].positive ? 1 : -1)*wl_best.evaluate(time_samples[j].time)) /
															(1 - mu*mu);
	}

	std::cout << "w=" << wl_best.my_weight << " ";
	wl_best.print();
	std::cout << std::endl;
	weak_classifiers.push_back(wl_best);
}

void Adaboost::save(FILE *file, bool lossy) {
	fwrite(&sum_alpha, sizeof(double), 1, file);
	int size = weak_classifiers.size();
	fwrite(&size, sizeof(int), 1, file);
	for (int i = 0; i < size; ++i) {
		weak_classifiers[i].save(file, lossy);
	}
}

void Adaboost::load(FILE *file) {
	fread(&sum_alpha, sizeof(double), 1, file);
	int size;
	fread(&size, sizeof(int), 1, file);
	weak_classifiers.resize(size);
	for (int i = 0; i < size; ++i) {
		weak_classifiers[i].load(file);
	}
}

void Adaboost::exportToArray(double *array, int maxLen, int &pos) {
	array[pos++] = sum_alpha;
	array[pos++] = weak_classifiers.size();
	for (int i = 0; i < weak_classifiers.size(); ++i) {
		weak_classifiers[i].exportToArray(array, maxLen, pos);
	}
}

void Adaboost::importFromArray(double *array, int len, int &pos) {
	sum_alpha = array[pos++];
	int size = array[pos++];
	weak_classifiers.resize(size);
	for (int i = 0; i < size; ++i) {
		weak_classifiers[i].importFromArray(array, len, pos);
	}
}

double Adaboost::classify(long time) const {
	double result = 0;
	for (int i = 0; i < weak_classifiers.size(); ++i) {
		result += weak_classifiers[i].evaluate(time) * weak_classifiers[i].my_weight;
	}
	return result/sum_alpha;
}

void Adaboost::print() {
	std::cout << "Adaboost [";
	for (int i = 0; i < weak_classifiers.size(); ++i) {
		weak_classifiers[i].print();
	}
	std::cout << "]" << std::endl;
}
