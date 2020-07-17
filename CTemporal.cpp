#include "CFrelement.h"
#include "CPerGaM.h"
#include "CTimeAdaptiveHist.h"
#include "CTimeHist.h"
#include "CTimeNone.h"
#include "CTimeMean.h"
//#include "CMises.h"
#include "CPythonHyperTime.h"
#include "CHyperTime.h"
#include "CExpectation.h"
#include "CMoments.h"
#include "CAdaboost.h"
#include "CNeural.h"
#include "CExpectationVM.h"
#include "CMomentsVM.h"
#include "CExpectationUW.h"

const char *temporalModelName[] = 
{
	"None",
	"Mean",
	"Hist",
	"FreMEn",
	"HyT-EM",
	"HyT-KM",
	"Gaussian",
	"Adaptive",
	"VonMises",
	"HyT-CEM-wC",
	"HyT-MM-wC",
	"Adaboost",
	"Neural",
	"HyT-CEM-vM",
	"HyT-MM-vM",
	"LiT-EM",
	"Number",
};


CTemporal* spawnTemporalModel(ETemporalType type,int maxPeriod,int elements,int numClasses)
{
	CTemporal *temporalModel;
	switch (type)
	{
		case TT_NONE: 		temporalModel = new CTimeNone(0);		break;
		case TT_MEAN: 		temporalModel = new CTimeMean(0);		break;
		case TT_HISTOGRAM:	temporalModel = new CTimeHist(0);		break;
		case TT_FREMEN: 	temporalModel = new CFrelement(0);		break;
		case TT_HYPER: 		temporalModel = new CHyperTime(0);		break;
		case TT_PYTHON: 	temporalModel = new CPythonHyperTime(0);	break;

		case TT_PERGAM: 	temporalModel = new CPerGaM(0);			break;
		case TT_ADAPTIVE: 	temporalModel = new CTimeAdaptiveHist(0);	break;
//		case TT_MISES: 		temporalModel = new CMises(0);			break;
		case TT_EXPECTATION:	temporalModel = new CExpectation(elements); break;
		case TT_MOMENTS:	temporalModel = new CMoments(elements);	break;
		case TT_ADABOOST:	temporalModel = new CAdaboost(elements); break;
		case TT_NEURAL:	temporalModel = new CNeural(elements); break;
		case TT_EXPECTATION_VM:	temporalModel = new CExpectationVM(elements); break;
		case TT_MOMENTS_VM:	temporalModel = new CMomentsVM(elements);	break;
		case TT_EXPECTATION_UW:	temporalModel = new CExpectationUW(elements); break;
		default: 		temporalModel = new CTimeNone(0);
	}
	temporalModel->init(maxPeriod,elements,numClasses);
	return temporalModel;
}

CTemporal* spawnTemporalModel(const char* type,int maxPeriod,int elements,int numClasses)
{
	int i = TT_NONE;
	for (i=0;i<TT_NUMBER && strcmp(type,temporalModelName[i])!=0;i++){}
	return spawnTemporalModel( (ETemporalType)i,maxPeriod,elements,numClasses);
}
