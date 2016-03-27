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
 * convert rgb color to equivalent shade using known weights
 * @param r - red portion of input color in range [0; 1]
 * @param g - green portion of input color in range [0; 1]
 * @param b - blue portion of input color in range [0; 1]
 * @return output shade in range [0; HISTOGRAM_LEVELS]
 */
double TMO2scaleToneManagement::RgbToGray(double r, double g, double b){
	return (0.299 * r + 0.587 * g + 0.114 * b) * HISTOGRAM_LEVELS;
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
void TMO2scaleToneManagement::FillHistogram(pfstmo::Array2D * image, int * histogram){
	int index = 0;	
	
	for (int j = 0; j < (*image).getRows(); j++){
		for (int i = 0; i < (*image).getCols(); i++){			
			index = (int) (*image)(i, j);
			
			// out of range fix
			if (index < 0) index = 0;
			if (index > (HISTOGRAM_LEVELS - 1)) index = HISTOGRAM_LEVELS - 1;
			
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
 * @param histogram - comulative histogram to look in
 * 
 * @return closest index 
 */
int TMO2scaleToneManagement::FindClosestShade(int value, int * histogram){		
	for (int i = 1; i < HISTOGRAM_LEVELS; i++){
		if (histogram[i] >= value){
			
			// find out if we are closer to current or previous item
			if ((histogram[i] - value) < (value - histogram[i-1])){
				//std::cerr << "FCS value " << value << ", index: " << i << std::endl;
				return i;
			}else{
				//std::cerr << "FCS value " << value << ", index: " << (i - 1) << std::endl;
				return (i - 1);
			}
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
	
	if (DEBUG) std::cerr << "NormaliseHistogram divider: " << divider << std::endl;
	
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
 * @param input - input image
 * @return weight
 */
double TMO2scaleToneManagement::BilateralFilterWeight(double i, double j, double k, double l, pfstmo::Array2D * input){		
	double numerator1 = pow(i - k, 2) + pow(j - l, 2);
	double denomirator1 = 2 * pow(sigmaS, 2);
	double numerator2 = pow((*input)(i, j) - (*input)(k, l), 2);	
	double denomirator2 = 2 * pow(sigmaR, 2);
	
	return exp(-(numerator1 / denomirator1) - (numerator2 / denomirator2));
}

/**
 * bilateral filter
 * 
 * @param output - resulted filtered image
 * @param input - input image
 */
void TMO2scaleToneManagement::BilateralFilter(pfstmo::Array2D * output, pfstmo::Array2D * input){
	double numerator = 0.0;			
	double denomirator = 0.0;
	double weight = 0.0;	
	
	// loop for each pixel in input image
	for (int j = 0; j < (*input).getRows(); j++){		
		for (int i = 0; i < (*input).getCols(); i++){
			numerator = 0.0;			
			denomirator = 0.0;
			weight = 0.0;
			
			// loop for each pixel in kernel
			for (int l = ((j - WINDOW_SIZE / 2) > 0) ? (j - WINDOW_SIZE / 2) : 0; l <= (j + WINDOW_SIZE / 2) && l < (*input).getRows(); l++){
				for (int k = ((i - WINDOW_SIZE / 2) > 0) ? (i - WINDOW_SIZE / 2) : 0; k <= (i + WINDOW_SIZE / 2) && k < (*input).getCols(); k++){
					weight = BilateralFilterWeight(i, j, k, l, input);					
					numerator += (*input)(k, l) * weight;
					denomirator += weight;					
				}
				
			}
			
			(*output)(i, j) = numerator / denomirator;
		}
	}		
}

/**
 * weight function for bilateral filter
 * 
 * @param i - x coordinate in input picture
 * @param j - y coordinate in input picture
 * @param k - x coordinate of pixel in neighbourhood
 * @param j - y coordinate of pixel in neighbourhood
 * @param input - input image
 * @param HPofInput - input filtered by highPass filter
 * @param numerator - if this is true then we are calculating numerator of BF, if this is false then we are calculating denumerator
 * @return weight
 */
double TMO2scaleToneManagement::CrossBilateralFilterWeight(double i, double j, double k, double l, pfstmo::Array2D * input, pfstmo::Array2D * HPofInput, bool numerator){	
	double sigmaSlocal = sigmaS * SIGMA_S_TEX_MULTIPLIER;
	double sigmaRlocal = sigmaR * SIGMA_S_TEX_MULTIPLIER;
	
	double numerator1 = pow(i - k, 2) + pow(j - l, 2);
	double denomirator1 = 2 * pow(sigmaSlocal, 2);
	//double numerator2 = (numerator) ? pow((*input)(i, j) - (*HPofInput)(k, l), 2) : pow((*input)(i, j) - (*input)(k, l), 2);	
	double numerator2 = (numerator) ? pow((*HPofInput)(i, j) - (*HPofInput)(k, l), 2) : pow((*input)(i, j) - (*input)(k, l), 2);	
	double denomirator2 = 2 * pow(sigmaRlocal, 2);
	
	return exp(-(numerator1 / denomirator1) - (numerator2 / denomirator2));
}

/**
 * @param output - filtered image
 * @param input - input image
 * @param HPofInput - input image filtered by highpass filter and passed throuhg abs()
 */
void TMO2scaleToneManagement::CrossBilateralFilter(pfstmo::Array2D * output, pfstmo::Array2D * input, pfstmo::Array2D * HPofInput){
	double numerator = 0.0;			
	double denomirator = 0.0;
	double weight = 0.0;	
	
	// loop for each pixel in input image
	for (int j = 0; j < (*input).getRows(); j++){		
		for (int i = 0; i < (*input).getCols(); i++){
			numerator = 0.0;			
			denomirator = 0.0;
			weight = 0.0;
			
			// loop for each pixel in kernel
			for (int l = ((j - WINDOW_SIZE / 2) > 0) ? (j - WINDOW_SIZE / 2) : 0; l <= (j + WINDOW_SIZE / 2) && l < (*input).getRows(); l++){
				for (int k = ((i - WINDOW_SIZE / 2) > 0) ? (i - WINDOW_SIZE / 2) : 0; k <= (i + WINDOW_SIZE / 2) && k < (*input).getCols(); k++){
					weight = CrossBilateralFilterWeight(i, j, k, l, input, HPofInput, true);
					numerator += (*HPofInput)(k, l) * weight;
					
					weight = CrossBilateralFilterWeight(i, j, k, l, input, HPofInput, false);
					denomirator += weight;					
				}
				
			}
			
			(*output)(i, j) = numerator / denomirator;
		}
	}			
}

/**
 * create grayscale image from RGB image based on color weights 
 * 
 * @param output - output grayscale image
 * @param input - input image
 */
void TMO2scaleToneManagement::CreateGrayscale(pfstmo::Array2D output, TMOImage * input){
	for (int j = 0; j < input->GetHeight(); j++){		
		for (int i = 0; i < input->GetWidth(); i++){
			output(i, j) = RgbToGray(input->GetPixel(i, j)[0], input->GetPixel(i, j)[1], input->GetPixel(i, j)[2]);			
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

/**
 * compute sigma S parameter for bilateral filter from imate width and height
 * 
 * @param width - width of picture in pixels
 * @param width - height of picture in pixels
 * @return computed sigmaS parameter
 */
double TMO2scaleToneManagement::ComputeSigmaS(int width, int height){
	return ((height < width) ? height : width) / 16.0;
}

/**
 * this method performs histogram matching
 * 
 * @param comulativeInputHistogram - normalised comulative histogram of input picutre (base)
 * @param comulativeModelHistogram - comulative histogram if model (model base) normalised to same value
 * @param input - array of values for histogram matching, this will beoverwritten by new values
 */
void TMO2scaleToneManagement::HistogramMatching(int * comulativeInputHistogram, int * comulativeModelHistogram, pfstmo::Array2D * input){
	int value, inputShade; 		
	
	for (int j = 0; j < (*input).getRows(); j++){
		for (int i = 0; i < (*input).getCols(); i++){
			inputShade = (int) (*input)(i, j);
			
			if (inputShade < 0) inputShade = 0;
			if (inputShade > (HISTOGRAM_LEVELS - 1)) inputShade = HISTOGRAM_LEVELS - 1;			
			
			value = comulativeInputHistogram[inputShade];
			(*input)(i, j) = FindClosestShade(value, comulativeModelHistogram);
			
			// DEBUG
			/*if ((i % 25 == 0) && (j % 25 == 0)){
				std::cerr << "IN HM inputShade: " << inputShade << std::endl;
				std::cerr << "IN HM value: " << value << std::endl;
				std::cerr << "IN HM result: " << input(i, j) << std::endl << std::endl;
			}*/
		}
	}
}

/**
 * fills the ro variable for textureness transfer
 * 
 * @param ro - ro field, output parameter
 * @param texturenessDesired - T' - input parameter
 * @param texturenessBase - T(B) - input parameter
 * @param texturenessDetail T(D) - input parameter
 */
void TMO2scaleToneManagement::FillRo(pfstmo::Array2D ro, pfstmo::Array2D texturenessDesired, pfstmo::Array2D texturenessBase, pfstmo::Array2D texturenessDetail){
	double max = std::numeric_limits<double>::min();
	
	for (int j = 0; j < texturenessDesired.getRows(); j++){
		for (int i = 0; i < texturenessDesired.getCols(); i++){			
			ro(i, j) = std::max(0.0f, (texturenessDesired(i, j) - texturenessBase(i, j)) / texturenessDetail(i, j));
			
			// find maximum
			if (ro(i, j) > max) max = ro(i, j);
		}
	}
	
	// normalise ro
	if (DEBUG) std::cerr << "FillRo max: " << max << std::endl;
	//double divider = (double) max / MAX_RHO;	
	//double divider = 10.0;
	
	//std::cerr << "FillRo divider: " << divider << std::endl;
	
	for (int j = 0; j < ro.getRows(); j++){
		for (int i = 0; i < ro.getCols(); i++){
			//ro(i, j) /= divider;
			
			// limit max value of rho
			if (ro(i, j) > MAX_RHO) ro(i, j) = MAX_RHO;
		}
	}
}

/**
 * high-pass filter
 * 
 * @param output - result: filtered array
 * @param input - input array
 * @param cutoff - cutoff
 */
void TMO2scaleToneManagement::HighPassFilter(pfstmo::Array2D output, pfstmo::Array2D input, double cutoff){		
	// VERSION 1 - ideal filter
	/*for (int j = 0; j < input.getRows(); j++){
		for (int i = 0; i < input.getCols(); i++){
			//output(i, j) = (input(i, j) > cutoff) ? input(i, j) : 0.0;
			output(i, j) = (input(i, j) > cutoff) ? input(i, j) : cutoff;
		}
	}*/	
	
	// VERSION 2 - Gaussian filter
	/*double RC = 1.0/(cutoff * 2 * M_PI);
	double alpha = RC / (RC + 1);	
	output(0, 0) = input(0, 0);
	
	// filtering
	for (int j = 0; j < input.getRows(); j++){
		for (int i = 0; i < input.getCols(); i++){
			int local_i = (i != 0) ? i : 1;
			int local_j = (i == 0 && j != 0) ? (j - 1) : j;
			
			output(local_i, local_j) = alpha * (output(local_i - 1, local_j) + input(local_i, local_j) - input(local_i - 1, local_j));
		}
	}*/
}

void TMO2scaleToneManagement::HighPassFilterV2(pfstmo::Array2D * input){		
	pfstmo::Array2D inputBackup = pfstmo::Array2D((*input).getCols(), (*input).getRows());
	
	// create backup array
	for (int j = 0; j < (*input).getRows(); j++){
		for (int i = 0; i < (*input).getCols(); i++){
			inputBackup(i, j) = (*input)(i, j);
		}
	}
	
	// VERSION 3 - 2D HP filter (http://homepages.inf.ed.ac.uk/rbf/BOOKS/PHILLIPS/cips2ed.pdf p.87)
	// filtering using this mask:
	// 0 -1 0
	// -1 5 -1
	// 0 -1 0
	for (int j = 1; j < (*input).getRows() - 1; j++){
		for (int i = 1; i < (*input).getCols() - 1; i++){
			/*(*input)(i, j) = 
				5.0 * inputBackup(i, j) + 		// current pixel
				(-1.0 * inputBackup(i - 1, j)) + 	// left pixel
				(-1.0 * inputBackup(i + 1, j)) + 	// right pixel
				(-1.0 * inputBackup(i, j + 1)) + 	// top pixel				
				(-1.0 * inputBackup(i, j - 1)); 	// bottom pixel*/
			
			(*input)(i, j) = inputBackup(i, j) * 5.0 - inputBackup(i + 1, j) - inputBackup(i - 1, j) - inputBackup(i, j - 1) - inputBackup(i, j + 1);

			// fix overflows
			if ((*input)(i, j) < 0.0) (*input)(i, j) = 0.0;
			if ((*input)(i, j) > MAX_VALUE_IN_RANGE) (*input)(i, j) = MAX_VALUE_IN_RANGE;
		}
	}
}

/**
 * computes textureness of given array using bilateral filter
 * assumes that sigmaS and sigmaR variables are set
 * 
 * @param texureness - output textureness
 * @param input - input array to compute textureness from
 * @param isDetail - indicate whether or not is this textureness of detail, detail has different scale and doesnt use high-pass filter
 */
void TMO2scaleToneManagement::FillTextureness(pfstmo::Array2D * textureness, pfstmo::Array2D * input, bool isDetail){		
	
	pfstmo::Array2D filteredInput = pfstmo::Array2D((*input).getCols(), (*input).getRows());
	//pfstmo::Array2D* filteredInput = new pfstmo::Array2D(input.getCols(),input.getRows());		
	
	// if (isDetail) -> simply copy input to filtered input
	for (int j = 0; j < (*input).getRows(); j++){
		for (int i = 0; i < (*input).getCols(); i++){				
			filteredInput(i, j) = (*input)(i, j) * 256.0;
		}
	}

	if (!isDetail){
		//HighPassFilter(filteredInput, input, sigmaS / SIGMA_S_TEX_MULTIPLIER);		
		HighPassFilterV2(&filteredInput);				
		
		// make absolute value of filtered array
		for (int j = 0; j < (*input).getRows(); j++){
			for (int i = 0; i < (*input).getCols(); i++){				
				filteredInput(i, j) = abs(filteredInput(i, j));
			}
		}
	}				
	
	CrossBilateralFilter(textureness, input, &filteredInput);		
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMO2scaleToneManagement::Transform(){
	// histograms for base
	int inputHistogram[HISTOGRAM_LEVELS];
	int modelHistogram[HISTOGRAM_LEVELS];
	int comulativeInputHistogram[HISTOGRAM_LEVELS];
	int comulativeModelHistogram[HISTOGRAM_LEVELS];	
	
	// histograms for textureness
	int inputHistogramTextureness[HISTOGRAM_LEVELS];
	int modelHistogramTextureness[HISTOGRAM_LEVELS];
	int comulativeInputHistogramTextureness[HISTOGRAM_LEVELS];
	int comulativeModelHistogramTextureness[HISTOGRAM_LEVELS];			
	
	// initialise histograms
	InitialiseHistogram(inputHistogram);
	InitialiseHistogram(comulativeInputHistogram);	
	InitialiseHistogram(modelHistogram);
	InitialiseHistogram(comulativeModelHistogram);
	InitialiseHistogram(inputHistogramTextureness);
	InitialiseHistogram(modelHistogramTextureness);	
	InitialiseHistogram(comulativeInputHistogramTextureness);
	InitialiseHistogram(comulativeModelHistogramTextureness);
	
	// get model filename
	std::string modelFilename = GetModelFilename(pSrc->GetFilename());
	if (modelFilename == ""){
		return -1;							// input filename doesnt contain _input. substring
	}

	// load model	
	model = new TMOImage(modelFilename.c_str());
	
	// bilateral filtering and tuxtureness variables
	pfstmo::Array2D base = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());			
	pfstmo::Array2D detail = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D grayscale = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D modelGrayscale = pfstmo::Array2D(model->GetWidth(),model->GetHeight());	
	pfstmo::Array2D modelBase = pfstmo::Array2D(model->GetWidth(),model->GetHeight());		
	pfstmo::Array2D texturenessInput = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D texturenessModel = pfstmo::Array2D(model->GetWidth(),model->GetHeight());
	pfstmo::Array2D texturenessDesired = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D texturenessBase = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D texturenessDetail = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D ro = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	
	// convert input image and model to rgb grayscale
	CreateGrayscale(grayscale, pSrc);
	CreateGrayscale(modelGrayscale, model);	
	
	// bilateralFiltering for input
	sigmaS = ComputeSigmaS(pSrc->GetWidth(), pSrc->GetHeight());
	sigmaR = 30.0;							// scale is the same is scale of pixel value [0, 255]		
	
	BilateralFilter(&base, &grayscale);
	GetDetailFromBase(detail, base, grayscale);				
	
	// bilateralFiltering for model	
	BilateralFilter(&modelBase, &modelGrayscale);	
	
	// get partial textureness	
	FillTextureness(&texturenessInput, &grayscale, false);	
	FillTextureness(&texturenessDesired, &grayscale, false);	
	FillTextureness(&texturenessModel, &modelGrayscale, false);	
	FillTextureness(&texturenessBase, &base, false);	
	FillTextureness(&texturenessDetail, &detail, true);
	
	// create histograms
	FillHistogram(&base, inputHistogram);	
	FillHistogram(&modelBase, modelHistogram);
	FillHistogram(&texturenessInput, inputHistogramTextureness);
	FillHistogram(&texturenessModel, modelHistogramTextureness);			
	
	// create comulative histograms
	ComputeComulativeHistogram(inputHistogram, comulativeInputHistogram);
	ComputeComulativeHistogram(modelHistogram, comulativeModelHistogram);
	ComputeComulativeHistogram(inputHistogramTextureness, comulativeInputHistogramTextureness);
	ComputeComulativeHistogram(modelHistogramTextureness, comulativeModelHistogramTextureness);
	
	// normalise comulative histograms
	NormaliseHistogram(comulativeInputHistogram, HISTOGRAM_NORMALISATION);	
	NormaliseHistogram(comulativeModelHistogram, HISTOGRAM_NORMALISATION);
	NormaliseHistogram(comulativeInputHistogramTextureness, HISTOGRAM_NORMALISATION);	
	NormaliseHistogram(comulativeModelHistogramTextureness, HISTOGRAM_NORMALISATION);	

	// DEBUG
	if (DEBUG){
		PrintHistogram(inputHistogram, "input base histogram");	
		PrintHistogram(modelHistogram, "model base histogram");	
		PrintHistogram(comulativeInputHistogram, "comulative input");
		PrintHistogram(comulativeModelHistogram, "comulative model");
		
		PrintHistogram(inputHistogramTextureness, "inputHistogramTextureness");	
		PrintHistogram(modelHistogramTextureness, "modelHistogramTextureness");	
		PrintHistogram(comulativeInputHistogramTextureness, "comulativeInputHistogramTextureness");
		PrintHistogram(comulativeModelHistogramTextureness, "comulativeModelHistogramTextureness");
	}
	
	// do historam matching from model base to new input base
	HistogramMatching(comulativeInputHistogram, comulativeModelHistogram, &base);
	
	// historam matching for textureness
	HistogramMatching(comulativeInputHistogramTextureness, comulativeModelHistogramTextureness, &texturenessDesired);		
	
	/*std::cerr << "DEBUG" << std::endl;
	std::cerr << "texturenessDesired(10, 10) " << texturenessDesired(10, 10) << std::endl;
	std::cerr << "texturenessDesired(50, 10) " << texturenessDesired(50, 10) << std::endl;
	std::cerr << "texturenessDesired(100, 100) " << texturenessDesired(100, 100) << std::endl;	
	std::cerr << "texturenessDesired(70, 70) " << texturenessDesired(70, 70) << std::endl;	
	std::cerr << "texturenessDesired(150, 150) " << texturenessDesired(150, 150) << std::endl;
	std::cerr << "texturenessDesired(210, 210) " << texturenessDesired(210, 210) << std::endl;	
	
	std::cerr << "texturenessInput(10, 10) " << texturenessInput(10, 10) << std::endl;
	std::cerr << "texturenessInput(50, 10) " << texturenessInput(50, 10) << std::endl;
	std::cerr << "texturenessInput(100, 100) " << texturenessInput(100, 100) << std::endl;	
	std::cerr << "texturenessInput(70, 70) " << texturenessInput(70, 70) << std::endl;	
	std::cerr << "texturenessInput(150, 150) " << texturenessInput(150, 150) << std::endl;
	std::cerr << "texturenessInput(210, 210) " << texturenessInput(210, 210) << std::endl;	
	
	std::cerr << "texturenessBase(10, 10) " << texturenessBase(10, 10) << std::endl;
	std::cerr << "texturenessBase(50, 10) " << texturenessBase(50, 10) << std::endl;
	std::cerr << "texturenessBase(100, 100) " << texturenessBase(100, 100) << std::endl;	
	std::cerr << "texturenessBase(70, 70) " << texturenessBase(70, 70) << std::endl;	
	std::cerr << "texturenessBase(150, 150) " << texturenessBase(150, 150) << std::endl;
	std::cerr << "texturenessBase(210, 210) " << texturenessBase(210, 210) << std::endl;	
	
	std::cerr << "texturenessDesired-texturenessBase(10, 10) " << texturenessDesired(10, 10)-texturenessBase(10, 10) << std::endl;
	std::cerr << "texturenessDesired-texturenessBase(50, 10) " << texturenessDesired(50, 10)-texturenessBase(50, 10) << std::endl;
	std::cerr << "texturenessDesired-texturenessBase(100, 100) " << texturenessDesired(100, 100)-texturenessBase(100, 100) << std::endl;	
	std::cerr << "texturenessDesired-texturenessBase(70, 70) " << texturenessDesired(70, 70)-texturenessBase(70, 70) << std::endl;	
	std::cerr << "texturenessDesired-texturenessBase(150, 150) " << texturenessDesired(150, 150)-texturenessBase(150, 150) << std::endl;
	std::cerr << "texturenessDesired-texturenessBase(210, 210) " << texturenessDesired(210, 210)-texturenessBase(210, 210) << std::endl;	
	
	std::cerr << "texturenessDetail(10, 10) " << texturenessDetail(10, 10) << std::endl;
	std::cerr << "texturenessDetail(50, 10) " << texturenessDetail(50, 10) << std::endl;
	std::cerr << "texturenessDetail(100, 100) " << texturenessDetail(100, 100) << std::endl;	
	std::cerr << "texturenessDetail(70, 70) " << texturenessDetail(70, 70) << std::endl;	
	std::cerr << "texturenessDetail(150, 150) " << texturenessDetail(150, 150) << std::endl;
	std::cerr << "texturenessDetail(210, 210) " << texturenessDetail(210, 210) << std::endl;*/
	
	// fill ro array
	FillRo(ro, texturenessDesired, texturenessBase, texturenessDetail);			
	
	// DEBUG
	/*InitialiseHistogram(inputHistogram);
	FillHistogram(base, inputHistogram);
	PrintHistogram(inputHistogram, "input base histogram after matching");	*/
	
	// show debug information
	if (DEBUG){
		std::cerr << "RO FOR DEBUG" << std::endl;
		for (int j = 0; j < ro.getRows(); j+=5){
			for (int i = 0; i < ro.getCols(); i+=5){
				//if (ro(i, j) > 100) std::cerr << "A ro is:" << ro(i, j) << std::endl;
				std::cerr << "ro(" << i << ", " << j << ") " << ro(i, j) << std::endl;
			}
		}
	}
	
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	double r, g, b, shade;
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
			
			// RGB setting
			
			// show base
			//shade = base(i, j) / HISTOGRAM_LEVELS;
			
			// show detail
			//shade = detail(i, j) / HISTOGRAM_LEVELS;
			
			// show base + detail
			//shade = (base(i, j) + detail(i, j) * 3) / HISTOGRAM_LEVELS;
						
			
			// show textureness
			//shade = texturenessInput(i, j) / MAX_VALUE_IN_RANGE;
			
			//shade = ro(i, j) / 5.0;
			
			// final result
			shade = (base(i, j) + ro(i, j) * detail(i, j)) / HISTOGRAM_LEVELS;			
			
			// fix "overflows"
			if (shade > 1.0){
				shade = 1.0;
			}else if (shade < 0.0){
				shade = 0.0;
			}
			
			*pDestinationData++ = shade;
			*pDestinationData++ = shade;
			*pDestinationData++ = shade;
		}
	}		
	
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);		
	
	return 0;
}

