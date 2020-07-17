#include <iostream>
#include <fstream>
#include "CMomentsVM.h"

using namespace std;

CMomentsVM::CMomentsVM(int idd) :
	CMoments(idd)
{
	type = TT_MOMENTS;

	for (int i = class_model.get_cluster_count(); i; --i) {
		densities.push_back(std::move(DensityParams::create(this, DensityParams::VON_MISES)));
		estimators.push_back(densities[densities.size()-1]->get_moment_estimator());
	}

	srandom(time(0));
}

CMomentsVM::~CMomentsVM()
{
}
