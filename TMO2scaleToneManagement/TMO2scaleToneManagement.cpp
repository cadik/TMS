/* --------------------------------------------------------------------------- *
 * TMO2scaleToneManagement.cpp: implementation of the TMO2scaleToneManagement class.   *
 * --------------------------------------------------------------------------- */

#include "TMO2scaleToneManagement.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMO2scaleToneManagement::TMO2scaleToneManagement()
{
	SetName(L"2scaleToneManagement");
	SetDescription(L"Method focused to pictorial looks, based on two-scale non-linear decomposition of an image. Input filename has to contain _input. substring, model filename has to contain _model substring.");
}

TMO2scaleToneManagement::~TMO2scaleToneManagement()
{
}

/**
 * initialise array with input histogram to 0 
 * 
 * @param histogram - array with histogram to initialise
 */
void TMO2scaleToneManagement::InitialiseHistogram(int * histogram){
	for (int i = 0; i < HISTOGRAM_LEVELS; i++){
		histogram[i] = 0;
	}
}

/**
 * print input histogram to stderr
 * 
 * @param histogram - array with histogram to fill
 * @param histogramName - name of histogram, will be printed first
 */
void TMO2scaleToneManagement::PrintHistogram(int * histogram, std::string histogramName){
	std::cerr << std::endl << "Histogram " << histogramName << ":" << std::endl;
	
	for (int i = 0; i < HISTOGRAM_LEVELS; i++){
		//std::cerr << "Histogram [" << i << "]: " << histogram[i] << std::endl;
		std::cerr << histogram[i] << std::endl;
	}
}

/**
 * fill the inputHistogram array with histogram data
 * 
 * @param image - input image to create histogram from
 * @param histogram - array for resulting histogram
 */
void TMO2scaleToneManagement::FillHistogram(TMOImage* image, int * histogram){
	int index = 0;	
	for (int j = 0; j < image->GetHeight(); j++){
		for (int i = 0; i < image->GetWidth(); i++){
			index = (0.299 * image->GetPixel(i, j)[0] + 0.587 * image->GetPixel(i, j)[1] + 0.114 * image->GetPixel(i, j)[2]) * HISTOGRAM_LEVELS;			
			histogram[index]++;
		}
	}
}

/**
 * finds filename of model based on filename of input
 * 
 * @param inputFilename - filename of input filename
 * @return filename of model
 */
std::string TMO2scaleToneManagement::GetModelFilename(std::string inputFilename){
	std::string inputSubstr = INPUT_FILENAME_SUBSTR;
	std::string modelSubstr = MODEL_FILENAME_SUBSTR;
	unsigned int inputPosition = inputFilename.find(inputSubstr);
	
	if (inputPosition == std::string::npos){
		std::cerr << "Error! Input filename has to contain _input. substring." << std::endl;
		return "";
	}
	
	return inputFilename.replace(inputPosition, inputSubstr.length(), modelSubstr);
}

/**
 * compute comulative histogram from classic histogram
 * 
 * @param histogram - input histogram
 * @param comulativeHistogram - output comulative histogram
 */
void TMO2scaleToneManagement::ComputeComulativeHistogram(int * histogram, int * comulativeHistogram){
	comulativeHistogram[0] = histogram[0];
	
	for (int i = 1; i < HISTOGRAM_LEVELS; i++){
		comulativeHistogram[i] = comulativeHistogram[i - 1] + histogram[i];
	}
}

/**
 * find on what index is given value in histogram
 * 
 * @param value - value to look for in histogram
 * @param histogram - histogram to look in
 * 
 * @return closest index 
 */
int TMO2scaleToneManagement::FindClosestShade(int value, int * histogram){
	for (int i = 1; i < HISTOGRAM_LEVELS; i++){
		if (histogram[i] >= value){
			return i;
		}
	}
	
	return HISTOGRAM_LEVELS - 1;
}

/**
 * normalise histogram (or comulative histogram)
 * 
 * @param histogram - histogram to normalise
 * @param value - value to normalise to
 */
void TMO2scaleToneManagement::NormaliseHistogram(int * histogram, int value){
	// find maximum
	int max = -1;
	for (int i = 0; i < HISTOGRAM_LEVELS; i++){
		if (histogram[i] > max) max = histogram[i];
	}
	
	// find divider
	double divider = (double) max / value;
	
	std::cerr << "divider: " << divider << std::endl;
	
	// normalise
	for (int i = 0; i < HISTOGRAM_LEVELS; i++){
		histogram[i] /= divider;
	}
}

/**
 * weight function for bilateral filter
 * 
 * @param i - x coordinate in input picture
 * @param j - y coordinate in input picture
 * @param k - x coordinate of pixel in neighbourhood
 * @param j - y coordinate of pixel in neighbourhood
 * @param sigmaS - filter parameters
 * @param sigmaR - filter parameters 
 * @param input - input image
 * @return weight
 */
double TMO2scaleToneManagement::BilateralFilterWeight(double i, double j, double k, double l, double sigmaS, double sigmaR, pfstmo::Array2D input){
	double numerator1 = pow(i - k, 2) + pow(j - l, 2);
	double denomirator1 = 2 * pow(sigmaS, 2);
	double numerator2 = abs(input(i, j) - input(k, l));
	double denomirator2 = 2 * pow(sigmaR, 2);
	
	return exp(-(numerator1 / denomirator1) - (numerator2 / denomirator2));
}

/**
 * bilateral filter
 * 
 * @param output - resulted filtered image
 * @param input - input image
 * @param sigmaS - filter parameters
 * @param sigmaR - filter parameters 
 */
void TMO2scaleToneManagement::BilateralFilter(pfstmo::Array2D output, pfstmo::Array2D input, double sigmaS, double sigmaR){	
	for (int j = 0; j < input.getRows(); j++){		
		for (int i = 0; i < input.getCols(); i++){
			double numerator = 0.0;			
			double denomirator = 0.0;
			double weight = 0.0;
						
			for (int l = ((j - WINDOW_SIZE / 2) > 0) ? (j - WINDOW_SIZE / 2) : 0; l < WINDOW_SIZE && l <= j; l++){
				for (int k = ((i - WINDOW_SIZE / 2) > 0) ? (i - WINDOW_SIZE / 2) : 0; k < WINDOW_SIZE && k <= i; k++){
					weight = BilateralFilterWeight(i, j, k, l, sigmaS, sigmaR, input);
					numerator += input(k, l) * weight;
					denomirator += weight;
				}
				
			}
			
			output(i, j) = numerator / denomirator;
		}
	}
}

/**
 * create grayscale image from RGB image based on color weights
 * for testing only?
 * 
 * @param output - output grayscale image
 * @param input - input image
 */
void TMO2scaleToneManagement::CreateGrayscale(pfstmo::Array2D output, TMOImage * input){
	for (int j = 0; j < input->GetHeight(); j++){		
		for (int i = 0; i < input->GetWidth(); i++){
			output(i, j) = 0.299 * input->GetPixel(i, j)[0] + 0.587 * input->GetPixel(i, j)[1] + 0.114 * input->GetPixel(i, j)[2];
		}
	}	
}

/**
 * simply deduct detail from base and original image
 * 
 * @param detail - output parameter, detail of the picture
 * @param base - input parameter, base of picture from bilateral filter
 * @param input - input image
 */
void TMO2scaleToneManagement::GetDetailFromBase(pfstmo::Array2D detail, pfstmo::Array2D base, pfstmo::Array2D input){
	for (int j = 0; j < input.getRows(); j++){		
		for (int i = 0; i < input.getCols(); i++){
			detail(i, j) = input(i, j) - base(i, j);
		}
	}		
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMO2scaleToneManagement::Transform(){
	int inputHistogram[HISTOGRAM_LEVELS];
	int modelHistogram[HISTOGRAM_LEVELS];
	int comulativeInputHistogram[HISTOGRAM_LEVELS];
	int comulativeModelHistogram[HISTOGRAM_LEVELS];
	
	// bilateral filtering variables
	pfstmo::Array2D base = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D detail = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D grayscale = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	
	// get model filename
	std::string modelFilename = GetModelFilename(pSrc->GetFilename());
	if (modelFilename == ""){
		return -1;							// input filename doesnt contain _input. substring
	}

	// load model
	// TODO remove TMOImage * ?
	TMOImage * model = new TMOImage(modelFilename.c_str());
	
	InitialiseHistogram(inputHistogram);
	InitialiseHistogram(comulativeInputHistogram);
	FillHistogram(pSrc, inputHistogram);	
	ComputeComulativeHistogram(inputHistogram, comulativeInputHistogram);
	NormaliseHistogram(comulativeInputHistogram, HISTOGRAM_NORMALISATION);	
	
	InitialiseHistogram(modelHistogram);
	InitialiseHistogram(comulativeModelHistogram);
	FillHistogram(model, modelHistogram);
	ComputeComulativeHistogram(modelHistogram, comulativeModelHistogram);
	NormaliseHistogram(comulativeModelHistogram, HISTOGRAM_NORMALISATION);

	PrintHistogram(inputHistogram, "input");	
	PrintHistogram(modelHistogram, "model");	
	PrintHistogram(comulativeInputHistogram, "comulative input");
	PrintHistogram(comulativeModelHistogram, "comulative model");	
	
	CreateGrayscale(grayscale, pSrc);
	
	// bilateralFiltering
	BilateralFilter(base, grayscale, 10, 10);
	GetDetailFromBase(detail, base, grayscale);		
	
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	double r, g, b;
	int inputShade;
	int value;
	int outputShade;
        
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++){
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++){
			r = *pSourceData++;
			g = *pSourceData++;
			b = *pSourceData++;
			
			// histogram matching
			inputShade = (0.299 * r + 0.587 * g + 0.114 * b) * HISTOGRAM_LEVELS;			
			value = comulativeInputHistogram[inputShade];			
			outputShade = FindClosestShade(value, comulativeModelHistogram);			
			
			// RGB setting
			*pDestinationData++ = (double) outputShade / HISTOGRAM_LEVELS;
			*pDestinationData++ = (double) outputShade / HISTOGRAM_LEVELS;
			*pDestinationData++ = (double) outputShade / HISTOGRAM_LEVELS;
			
			// LAB setting
			/**pDestinationData++ = base(i, j);
			*pDestinationData++ = 0.0;
			*pDestinationData++ = 0.0;*/
		}
	}
	
	pSrc->ProgressBar(j, pSrc->GetHeight());
	
	pDst->Convert(TMO_RGB);
	return 0;
}

