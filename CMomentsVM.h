#ifndef CMomentsVM_H
#define CMomentsVM_H

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <memory>
#include "CTemporal.h"
#include "CMoments.h"

#include "moments/moment_estimator.h"
#include "moments/density_params.h"

#define MAX_ID_LENGTH 100
	
using namespace std;

class CMomentsVM: public CMoments
{
	public:
		CMomentsVM(int idd);
		~CMomentsVM();
};

#endif
