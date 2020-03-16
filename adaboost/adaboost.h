#ifndef ADABOOST_H
#define ADABOOST_H

#include <vector>
#include "logreg.hpp"

class Adaboost
{
	public:
		Adaboost();
		void add_time(long int time, bool positive);
		void add_period(double period);
		void train();
		double classify(long int time) const;

		void save(FILE* file, bool lossy = false);
		void load(FILE* file);
		void exportToArray(double* array, int maxLen, int& pos);
		void importFromArray(double* array, int len, int& pos);

		void print();

	private:

		class TimeSample {
			public:
				TimeSample() :
					time(),
					weight(),
					positive()
				{ }
				TimeSample(long int time_, float weight_, bool positive_) :
					time(time_),
					weight(weight_),
					positive(positive_)
				{ }
				double time;
				double weight;
				bool positive;
		};

		std::vector<TimeSample> time_samples;
		std::vector<double> periods;
		std::vector<Logreg> weak_classifiers;
		long int min_time;
		long int max_time;
		double sum_alpha;

		void adaboost_iteration();
};

#endif // ADABOOST_H
