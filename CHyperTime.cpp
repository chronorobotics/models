#include "CHyperTime.h"

using namespace std;

CHyperTime::CHyperTime(int id)
{
	type = TT_HYPER;
	delete modelPositive;
	delete modelNegative;
//	modelPositive.release();
//	modelNegative.release();
	spaceDimension = 1;
	timeDimension = 0;
	maxTimeDimension = 10;
	covarianceType = EM::COV_MAT_GENERIC;
	positives = negatives = 0;
	corrective = 1.0;
	integral = 0;
}

void CHyperTime::init(int iMaxPeriod, int elements, int numClasses)
{
	maxPeriod = iMaxPeriod;
	numElements = elements;
}

CHyperTime::~CHyperTime()
{
}

// adds new state observations at given times
int CHyperTime::add(uint32_t time, float state)
{
	sampleArray[numSamples].t = time;
	sampleArray[numSamples].v = state;
	integral +=state;
	numSamples++;
	return 0; 
}

void CHyperTime::reinit_models_if_null() {
	if (modelPositive.empty()) {
		modelPositive = EM::create();
		modelPositive->setClustersNumber(order);
		modelPositive->setCovarianceMatrixType(covarianceType);
	}
	if (modelNegative.empty()) {
		modelNegative = EM::create();
		modelNegative->setClustersNumber(order);
		modelNegative->setCovarianceMatrixType(covarianceType);
	}
}

/*required in incremental version*/
void CHyperTime::update(int modelOrder,unsigned int* times,float* signal,int length)
{
	int numTraining = numSamples;
	int numEvaluation = 0;
	if (order != modelOrder) {
		modelPositive.release();
		modelNegative.release();
		order = modelOrder;
	}
	reinit_models_if_null();

	/*separate positives and negative examples*/
	Mat samplesPositive(positives, spaceDimension + timeDimension, CV_32FC1);
	Mat samplesNegative(negatives, spaceDimension + timeDimension, CV_32FC1);
	
	float vDummy = 0.5;
	long int tDummy = 0.5;
	for (int i = 0; i < numTraining; i++){
		vDummy = sampleArray[i].v;
		samplesPositive.push_back(vDummy);
		positives++;
		samplesNegative.push_back(vDummy);
		negatives++;
	}

	periods.clear();
	bool stop = false;
	do { 
		/*find the gaussian mixtures*/
		if (positives <= order || negatives <= order) break;
		modelPositive->trainEM(samplesPositive);
		modelNegative->trainEM(samplesNegative);
		if (modelNegative->isTrained()) std::cout << " Positive model is trained" << std::endl;
		if (modelPositive->isTrained()) std::cout << " Negative model is trained" << std::endl;
		std::cout << "Model trained with "<< order <<" clusters, "<< timeDimension <<"dimensions, "<< positives <<" positives and "<< negatives <<" negatives" << std::endl;
		print();

		/*analyse model error for periodicities*/
		CFrelement fremen(0);
		float err = 0;
		float sumErr = 0;
		fremen.init(maxPeriod, maxTimeDimension, 1);

		/*calculate model error across time*/
		for (int i = 0; i < numTraining; i++)
		{
			fremen.add(sampleArray[i].t, estimate(sampleArray[i].t)-sampleArray[i].v);
		}

		/*determine model weights*/
		float integralMod = 0;
		numEvaluation = numSamples;
		for (int i = 0; i < numEvaluation;i++) integralMod += estimate(sampleArray[i].t);
		corrective = corrective*integral/integralMod;
		
		/*calculate evaluation error*/
		for (int i = 0; i < numEvaluation; i++)
		{
			err = estimate(sampleArray[i].t)-sampleArray[i].v;
			sumErr+=err*err;
		}
		sumErr=sqrt(sumErr/numEvaluation);

		/*retrieve dominant error period*/	
		int maxOrder = 1;
		fremen.update(timeDimension/2+1);
		int period = fremen.getPredictFrelements()[0].period;
		bool expand = true;
		fremen.print(true);
		std::cout << "Model error with "<< timeDimension <<" time dimensions and "<< order <<" clusters is "<< sumErr << std::endl;

		/*if the period already exists, then skip it*/
		for (int d = 0; d < timeDimension/2; d++)
		{
			if (period == periods[d]) period = fremen.getPredictFrelements()[d+1].period;
		}
		errors[timeDimension/2] = sumErr;

		/*error has increased: cleanup and stop*/
		if (timeDimension > 1 && errors[timeDimension/2-1] < sumErr)
		{
			printf("Error increased from %.3f to %.3f\n", errors[timeDimension/2-1], errors[timeDimension/2]);
			timeDimension-=2;
			load("model");
			samplesPositive = samplesPositive.colRange(0, samplesPositive.cols-2);
			samplesNegative = samplesNegative.colRange(0, samplesNegative.cols-2);
			if (order < maxOrder){
				modelPositive.release();
				modelNegative.release();
			}
			reinit_models_if_null();
			std::cout << "Reducing hypertime dimension to "<< timeDimension <<": ";
			for (int i = 0; i < timeDimension/2; i++) std::cout << " " << periods[i];
			std::cout << std::endl;
			stop = true;
		}else{
			save("model");
		}
		if (timeDimension >= maxTimeDimension) stop = true;

		/*hypertime expansion*/
		if (stop == false && expand == true){
			printf("Adding period %i \n",period);
			Mat hypertimePositive(positives, 2, CV_32FC1);
			Mat hypertimeNegative(negatives, 2, CV_32FC1);
			positives = negatives = 0;
			for (int i = 0;i<numTraining;i++)
			{
				vDummy = sampleArray[i].v;
				tDummy = sampleArray[i].t;
				hypertimePositive.at<float>(positives,0)=cos((float)tDummy/period*2*M_PI);
				hypertimePositive.at<float>(positives,1)=sin((float)tDummy/period*2*M_PI);
				positives++;
				hypertimeNegative.at<float>(negatives,0)=cos((float)tDummy/period*2*M_PI);
				hypertimeNegative.at<float>(negatives,1)=sin((float)tDummy/period*2*M_PI);
				negatives++;
			}
			hconcat(samplesPositive, hypertimePositive, samplesPositive);
			hconcat(samplesNegative, hypertimeNegative, samplesNegative);
			periods.push_back(period);
			timeDimension += 2;
		}
		if (order <  maxOrder) stop = false;
	} while (!stop);
}

float CHyperTime::estimate(uint32_t t)
{
	/*is the model valid?*/
	if (modelNegative->isTrained() && modelPositive->isTrained()){
		Mat sample(1, spaceDimension + timeDimension, CV_32FC1);
		sample.at<float>(0,0) = 1;

		/*augment data sample with hypertime dimensions)*/
		for (int i = 0; i < timeDimension/2; i++){
			sample.at<float>(0, spaceDimension + 2*i + 0) = cos((float)t/periods[i]*2*M_PI);
			sample.at<float>(0, spaceDimension + 2*i + 1) = sin((float)t/periods[i]*2*M_PI);
		}
		Mat probs(1, 2, CV_32FC1);
		Vec2f a = modelPositive->predict2(sample, probs);

		sample.at<float>(0,0) = 0;
		Vec2f b = modelNegative->predict2(sample, probs);

		double d = ((positives*exp(a(0))+negatives*exp(b(0))));
		if (d > 0) return corrective*positives*exp(a(0))/d;
	} else {
		std::cout << "Model estimation skipped" << std::endl;
	}
	/*any data available?*/
	if (negatives+positives > 0) return (float)positives/(positives+negatives);
	return 0.5;
}

void CHyperTime::print(bool all)
{
	Mat meansPositive = modelPositive->getMeans();
	Mat meansNegative = modelNegative->getMeans();
	std::cout << meansPositive << std::endl;
	std::cout << meansNegative << std::endl;
	//std::cout << periods << std::endl;	
	std::cout << order <<" "<< timeDimension << std::endl;
}

float CHyperTime::predict(uint32_t time)
{
	return estimate(time);	
}

int CHyperTime::save(const char* name,bool lossy)
{
	FileStorage fsp(name, FileStorage::WRITE);
	fsp << "periods" << periods;
	fsp << "order" << order;
	fsp << "positives" << positives;
	fsp << "negatives" << negatives;
	fsp << "corrective" << corrective;
	cvStartWriteStruct(*fsp, "ModelPositive", CV_NODE_MAP);
	if (modelPositive->isTrained()){
		modelPositive->write(fsp);
		printf("saving positive\n");
	}
	cvEndWriteStruct(*fsp);
	cvStartWriteStruct(*fsp, "ModelNegative", CV_NODE_MAP);
	if (modelNegative->isTrained()){
		modelNegative->write(fsp);
		printf("saving negative\n");
	}
	cvEndWriteStruct(*fsp);
	fsp.release();

	return 0;
}

int CHyperTime::load(const char* name)
{
	FileStorage fs(name, FileStorage::READ);
	fs["periods"] >> periods;
	fs["order"] >> order;
	fs["positives"] >> positives;
	fs["negatives"] >> negatives;
	fs["corrective"] >> corrective;

	modelPositive = EM::create();
	modelNegative = EM::create();
	modelPositive->setClustersNumber(order);
	modelNegative->setClustersNumber(order);
	modelPositive->setCovarianceMatrixType(covarianceType);
	modelNegative->setCovarianceMatrixType(covarianceType);

	modelPositive->read(fs["ModelPositive"]);
	modelNegative->read(fs["ModelNegative"]);
	fs.release();

	timeDimension = periods.size()*2;

	return 0;
}

int CHyperTime::save(FILE* file,bool lossy)
{
//	int frk = numElements;
//	fwrite(&frk,sizeof(uint32_t),1,file);
//	fwrite(&storedGain,sizeof(float),1,file);
//	fwrite(storedFrelements,sizeof(SFrelement),numElements,file);
	return 0;
}

int CHyperTime::load(FILE* file)
{
//	int frk = numElements;
//	fwrite(&frk,sizeof(uint32_t),1,file);
//	fwrite(&storedGain,sizeof(float),1,file);
//	fwrite(storedFrelements,sizeof(SFrelement),numElements,file);
	return 0;
}

/*this is very DIRTY, but I don't see any other way*/
int CHyperTime::exportToArray(double* array,int maxLen)
{
	memset(array,0,sizeof(double)*maxLen);
	array[0] = TT_HYPER;
	if (modelNegative->isTrained() && modelPositive->isTrained()){
		save("hypertime.tmp");
		FILE*  file = fopen("hypertime.tmp","r");
		int len = fread(&array[3],1,maxLen,file);
		fclose(file);
		array[1] = len/sizeof(double)+1;

		return array[1]+3;
	}else{
		array[1] = 0;
		array[2] = positives;
		array[3] = negatives;
		array[4] = order;
		return 5;
	}
}

/*this is very DIRTY, but I don't see any other way*/
int CHyperTime::importFromArray(double* array,int len)
{
	if (array[1] > 0){
		FILE*  file = fopen("hypertime.tmp","w");
		fwrite(&array[3],1,(array[1])*sizeof(double),file);
		fclose(file);
		load("hypertime.tmp");
	}else{
		periods.clear();
		positives = array[2];
		negatives = array[3];
		order = array[4];
		modelPositive.release();
		modelNegative.release();
		reinit_models_if_null();
	}
	return 0;
}
