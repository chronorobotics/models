#include "CTimeHist.h"

using namespace std;

CTimeHist::CTimeHist(int idd)
{
	id=idd;
	firstTime = -1;
	lastTime = -1;
	measurements = 0;
	maxPeriod = 0;
	numElements = 0;
	correction = 0;
	type = TT_HISTOGRAM;
}

void CTimeHist::init(int iMaxPeriod,int elements,int models)
{
	maxPeriod = iMaxPeriod;
	maxPeriod = 86400;
	numElements = elements;
	numModels = models;
	predictHistogram = new float[numElements];
	storedHistogram = new float[numElements];
	measurementsHistogram = new int[numElements];
	for (int i=0;i<numElements;i++){
		predictHistogram[i] = 0.5; 
		measurementsHistogram[i] = 0;
		storedHistogram[i] = 0;
	}
}

CTimeHist::~CTimeHist()
{
	delete[] predictHistogram;
	delete[] measurementsHistogram;
	delete[] storedHistogram;
}

// adds new state observations at given times
int CTimeHist::add(uint32_t time,float state)
{
	if (measurements == 0) firstTime = time;
	lastTime = time;
	int index = ((time%maxPeriod)*numElements/maxPeriod)%numElements;
	storedHistogram[index] += state;
	measurementsHistogram[index] ++;
	measurements++;

	return 0; 
}

/*not required in incremental version*/
void CTimeHist::update(int modelOrder,unsigned int* times,float* signal,int length)
{
	for (int i=0;i<numElements;i++){
		 if (measurementsHistogram[i] > 0) predictHistogram[i] = storedHistogram[i]/measurementsHistogram[i]; else predictHistogram[i] = 0.5;
	}
}

/*text representation of the fremen model*/
void CTimeHist::print(bool verbose)
{
	std::cout << "Model " << id << " Size: " << measurements << " " << " Bins: " << numElements << " ";
	if (verbose){
		std::cout << "Bin values: ";
		for (int i = 0;i<numElements;i++) std::cout << storedHistogram[i] << " ";
	}
	std::cout << std::endl;
}

float CTimeHist::estimate(uint32_t time)
{
	int index = ((time%maxPeriod)*numElements/maxPeriod)%numElements;
	float estimate = 1.0/numModels;
	if (measurementsHistogram[index] > 0) estimate = storedHistogram[index]/measurementsHistogram[index];
	float saturation = 0.001;
	if (estimate > 1.0-saturation) estimate =  1.0-saturation;
	if (estimate < 0.0+saturation) estimate =  0.0+saturation;
	return estimate;
}

float CTimeHist::predict(uint32_t time)
{
	float estimate = predictHistogram[((time%maxPeriod)*numElements/maxPeriod)%numElements];
	float saturation = 0.001;
	if (estimate > 1.0-saturation) estimate =  1.0-saturation;
	if (estimate < 0.0+saturation) estimate =  0.0+saturation;
	return estimate;
}
int CTimeHist::save(const char* name,bool lossy)
{
	FILE* file = fopen(name,"w");
	save(file);
	fclose(file);
	return 0;
}

int CTimeHist::load(const char* name)
{
	FILE* file = fopen(name,"r");
	load(file);
	fclose(file);
	return 0;
}


int CTimeHist::save(FILE* file,bool lossy)
{
	double *array = new double[MAX_TEMPORAL_MODEL_SIZE*sizeof(double)];
	int len = exportToArray(array,MAX_TEMPORAL_MODEL_SIZE);
	fwrite(array,sizeof(double),len,file);
	delete[] array;
	return -1;
}

int CTimeHist::load(FILE* file)
{
	double *array = new double[MAX_TEMPORAL_MODEL_SIZE*sizeof(double)];
	int len = fread(array,sizeof(double),MAX_TEMPORAL_MODEL_SIZE,file);
	importFromArray(array,len);
	delete[] array;
	return len;
}

int CTimeHist::exportToArray(double* array,int maxLen)
{
	int pos = 0;
	array[pos++] = type;
	array[pos++] = numElements;
	array[pos++] = numModels;
	array[pos++] = id;
	array[pos++] = measurements;
	for (int i = 0;i<numElements && pos < MAX_TEMPORAL_MODEL_SIZE;i++) array[pos++] = (double)storedHistogram[i];
	for (int i = 0;i<numElements && pos < MAX_TEMPORAL_MODEL_SIZE;i++) array[pos++] = (double)measurementsHistogram[i];

	if (pos == MAX_TEMPORAL_MODEL_SIZE) std::cout << "Could not save the model dur to its size" << std::endl;
	return pos;
}

int CTimeHist::importFromArray(double* array,int len)
{
	int pos = 0;
	type = (ETemporalType)array[pos++];
	if (type != TT_HISTOGRAM) std::cerr << "Error loading the model, type mismatch." << std::endl;
	numElements = array[pos++];
	numModels = array[pos++];
	id = array[pos++];  
	measurements = array[pos++]; 
 
	delete[] storedHistogram;
	delete[] measurementsHistogram;
	init(maxPeriod,numElements,numModels);

	for (int i = 0;i<numElements && pos < MAX_TEMPORAL_MODEL_SIZE;i++)storedHistogram[i]=array[pos++]; 
	for (int i = 0;i<numElements && pos < MAX_TEMPORAL_MODEL_SIZE;i++)measurementsHistogram[i]=array[pos++]; 
	update(numElements);
	if (pos == MAX_TEMPORAL_MODEL_SIZE) std::cout << "Model was not properly saved before." << std::endl;
	return pos;
}
