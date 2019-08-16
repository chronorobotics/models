#ifndef EM_MULTI_ALTERNATIVE_H
#define EM_MULTI_ALTERNATIVE_H

#include <vector>

class EMMultiAlternative
{
	public:
		EMMultiAlternative(int cluster_count_, int dimension_);

		void train();
		void add_value(std::vector<bool> value);

		void print();

		std::vector<double> get_alpha_at(std::vector<bool> value) const;
		std::vector<std::vector<double> > get_means() const;
		int get_cluster_count() const;

	private:
		void expectation();
		double maximisation();

		int cluster_count;
		int dimension;

		class Cluster {
			public:
				Cluster(int dimension_);
				Cluster(std::vector<double> mean_, double weight_);
				std::vector<double> mean;
				double weight;

				void print();
				double probability_at(std::vector<bool> value) const;
		};

		class Measurement {
			public:
				Measurement(std::vector<bool> value_, int cluster_count);
				std::vector<bool> value;
				std::vector<double> alpha;
		};

		std::vector<Cluster> clusters;
		std::vector<Measurement> measurements;
};

#endif // EM_MULTI_ALTERNATIVE_H
