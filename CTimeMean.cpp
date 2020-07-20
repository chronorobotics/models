#include "CTimeMean.h"

using namespace std;

CTimeMean::CTimeMean(int idd)
{
	id=idd;
	firstTime = -1;
	lastTime = -1;
	measurements = 0;
	maxPeriod = 0;
	numElements = 0;
	positive = 0;
	correction = 0;
	type = TT_MEAN;
}

void CTimeMean::init(int iMaxPeriod,int elements,int numClasses)
{
	maxPeriod = iMaxPeriod;
	numElements = 1;
	estimation = 1.0/numClasses; 
}

CTimeMean::~CTimeMean()
{
}


// adds new state observations at given times
int CTimeMean::add(uint32_t time,float state)
{
	if (measurements == 0) firstTime = time;
	lastTime = time;
	positive=positive+state;
	measurements++;
	return 0; 
}

/*not required in incremental version*/
void CTimeMean::update(int modelOrder,unsigned int* times,float* signal,int length)
{
	if (measurements  > 0) estimation = positive/measurements;
}

/*text representation of the fremen model*/
void CTimeMean::print(bool verbose)
{
	std::cout << "Model " << id << " Size: " << measurements << " ";
	if (verbose){
		std::cout << "Mean: ";
		std::cout << (positive/measurements);
	}
	std::cout << std::endl;
}

float CTimeMean::estimate(uint32_t time)
{
	float estimate = estimation;//positive/measurements; 
	float saturation = 0.001;
	if (estimate > 1.0-saturation) estimate =  1.0-saturation;
	if (estimate < 0.0+saturation) estimate =  0.0+saturation;
	return estimate;
}

float CTimeMean::predict(uint32_t time)
{
	float estimate = estimation; 
	float saturation = 0.001;
	if (estimate > 1.0-saturation) estimate =  1.0-saturation;
	if (estimate < 0.0+saturation) estimate =  0.0+saturation;
	return estimate;
}
int CTimeMean::save(const char* name,bool lossy)
{
	FILE* file = fopen(name,"w");
	save(file);
	fclose(file);
	return 0;
}

int CTimeMean::load(const char* name)
{
	FILE* file = fopen(name,"r");
	load(file);
	fclose(file);
	return 0;
}


int CTimeMean::save(FILE* file,bool lossy)
{
	double array[10000];
	int len = exportToArray(array,10000);
	fwrite(array,sizeof(double),len,file);
	return 0;
}

int CTimeMean::load(FILE* file)
{
	double* array = new double[MAX_TEMPORAL_MODEL_SIZE];
	int len = fread(array,sizeof(double),MAX_TEMPORAL_MODEL_SIZE,file);
	importFromArray(array,len);
	delete[] array;
	return 0;
}

int CTimeMean::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = type;
	array[pos++] = positive;
	array[pos++] = id;
	array[pos++] = measurements;
	return pos;
}

int CTimeMean::importFromArray(double* array,int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_MEAN) std::cerr << "Error loading the model, type mismatch." << std::endl;
	positive = array[pos++];
	id = array[pos++];  
	measurements = array[pos++];
	update(0);
	return pos;
}
