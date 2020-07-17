#include <iostream>
#include <fstream>
#include "CMomentsVM.h"

using namespace std;

CMomentsVM::CMomentsVM(int idd) :
	CMoments(idd)
{
	type = TT_MOMENTS_VM;

  pos_density = DensityParams::create(this, DensityParams::SQDIST);
  neg_density = DensityParams::create(this, DensityParams::SQDIST);

	pos_estimator = pos_density->get_moment_estimator();
	neg_estimator = neg_density->get_moment_estimator();
}

CMomentsVM::~CMomentsVM()
{
}
