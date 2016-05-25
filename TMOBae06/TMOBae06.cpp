/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*                                 diploma thesis                               *
*             Author: Petr Pospisil [xpospi68 AT stud.fit.vutbr.cz]            *
*                                    Brno 2016                                 *
*                                                                              *
*******************************************************************************/

/* --------------------------------------------------------------------------- *
 * TMOBae06.cpp: implementation of the TMOBae06 class.                         *
 *               Two-scale Tone Management for Photographic Look               *
 * Method number: 3                                                            *
 * --------------------------------------------------------------------------- */

#include "TMOBae06.h"

/**
 * constructor, prepare parameters
 */
TMOBae06::TMOBae06()
{
	SetName(L"Bae06");
	SetDescription(L"Two-scale Tone Management for Photographic Look");
	
	// show/hide debug info
	verbose.SetName(L"v");
	verbose.SetDescription(L"Verbose output");
	verbose.SetDefault(false);
	verbose=false;	
	this->Register(verbose);
	
	// turn on special hdr settings
	hdr.SetName(L"HDR");
	hdr.SetDescription(L"enable/disable special HDR setting");
	hdr.SetDefault(false);
	hdr=false;	
	this->Register(hdr);	
	
	// model filename
	modelFileNameParam.SetName(L"model");
	modelFileNameParam.SetDescription(L"filename of model file - mandatory parameter");
	modelFileNameParam.SetDefault("");
	modelFileNameParam="";
	this->Register(modelFileNameParam);
}

TMOBae06::~TMOBae06()
{
}

/**
 * convert rgb color to equivalent shade using known weights
 * 
 * not used
 * 
 * @param r - red portion of input color in range [0; 1]
 * @param g - green portion of input color in range [0; 1]
 * @param b - blue portion of input color in range [0; 1]
 * @return output shade in range [0; HISTOGRAM_LEVELS]
 */
double TMOBae06::RgbToGray(double r, double g, double b){
	return (0.299 * r + 0.587 * g + 0.114 * b) * HISTOGRAM_LEVELS_F;
}

/**
 * initialise array with histogram to 0 
 * 
 * @param histogram - array with histogram to initialise
 */
void TMOBae06::InitialiseHistogram(int * histogram){
	for (int i = 0; i < HISTOGRAM_LEVELS; i++){
		histogram[i] = 0;
	}
}

/**
 * print histogram to stderr
 * 
 * @param histogram - array with histogram to fill
 * @param histogramName - name of histogram, will be printed first
 */
void TMOBae06::PrintHistogram(int * histogram, std::string histogramName){
	if (verbose) std::cerr << std::endl << "Histogram " << histogramName << ":" << std::endl;
	
	for (int i = 0; i < HISTOGRAM_LEVELS; i++){
		//if (verbose) std::cerr << "Histogram [" << i << "]: " << histogram[i] << std::endl;
		if (verbose) std::cerr << histogram[i] << std::endl;
	}
}

/**
 * fill the inputHistogram array with histogram data
 * 
 * @param image - input image to create histogram from
 * @param histogram - array for resulting histogram
 */
void TMOBae06::FillHistogram(pfstmo::Array2D * image, int * histogram){
	int index = 0;	
	
	for (int j = 0; j < (*image).getRows(); j++){
		for (int i = 0; i < (*image).getCols(); i++){			
			index = (int) ((*image)(i, j) * HISTOGRAM_BOOSTER_F);
			
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
 * removed
 * 
 * @param inputFilename - filename of input filename
 * @return filename of model
 */
/*std::string TMOBae06::GetModelFilename(std::string inputFilename){
	std::string inputSubstr = INPUT_FILENAME_SUBSTR;
	std::string modelSubstr = MODEL_FILENAME_SUBSTR;
	unsigned int inputPosition = inputFilename.find(inputSubstr);
	
	if (inputPosition == std::string::npos){
		if (verbose) std::cerr << "Error! Input filename has to contain _input. substring." << std::endl;
		return "";
	}
	
	return inputFilename.replace(inputPosition, inputSubstr.length(), modelSubstr);
}*/

/**
 * compute comulative histogram from classic histogram
 * 
 * @param histogram - input histogram
 * @param comulativeHistogram - output comulative histogram
 */
void TMOBae06::ComputeComulativeHistogram(int * histogram, int * comulativeHistogram){
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
int TMOBae06::FindClosestShade(int value, int * histogram){		
	for (int i = 1; i < HISTOGRAM_LEVELS; i++){
		if (histogram[i] >= value){
			
			// find out if we are closer to current or previous item
			if ((histogram[i] - value) < (value - histogram[i-1])){				
				return i;
			}else{
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
void TMOBae06::NormaliseHistogram(int * histogram, int value){
	// find maximum
	int max = -1;
	for (int i = 0; i < HISTOGRAM_LEVELS; i++){
		if (histogram[i] > max) max = histogram[i];
	}
	
	// find divider
	double divider = (double) max / value;
	
	if (verbose) std::cerr << "NormaliseHistogram divider: " << divider << std::endl;
	
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
double TMOBae06::BilateralFilterWeight(double i, double j, double k, double l, pfstmo::Array2D * input){		
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
void TMOBae06::BilateralFilter(pfstmo::Array2D * output, pfstmo::Array2D * input){
	double numerator = 0.0;			
	double denomirator = 0.0;
	double weight = 0.0;
	int windowSize = (hdr) ? WINDOW_SIZE_HDR : WINDOW_SIZE_LDR;
	
	// loop for each pixel in input image
	for (int j = 0; j < (*input).getRows(); j++){		
		for (int i = 0; i < (*input).getCols(); i++){
			numerator = 0.0;			
			denomirator = 0.0;
			weight = 0.0;
			
			// loop for each pixel in kernel
			for (int l = ((j - windowSize / 2) > 0) ? (j - windowSize / 2) : 0; l <= (j + windowSize / 2) && l < (*input).getRows(); l++){
				for (int k = ((i - windowSize / 2) > 0) ? (i - windowSize / 2) : 0; k <= (i + windowSize / 2) && k < (*input).getCols(); k++){
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
 * @param numerator - if this is true then we are calculating numerator of CBF, if this is false then we are calculating denumerator
 * @return weight
 */
double TMOBae06::CrossBilateralFilterWeight(double i, double j, double k, double l, pfstmo::Array2D * input, pfstmo::Array2D * HPofInput, bool numerator){	
	double sigmaSlocal = sigmaS * SIGMA_S_TEX_MULTIPLIER;
	double sigmaRlocal = sigmaR;
	
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
void TMOBae06::CrossBilateralFilter(pfstmo::Array2D * output, pfstmo::Array2D * input, pfstmo::Array2D * HPofInput){
	double numerator = 0.0;			
	double denomirator = 0.0;
	double weight = 0.0;
	
	// localWindowSize has to be bigger, but too big would take very long time
	int localWindowSize = (hdr) ? (WINDOW_SIZE_HDR * 2) : (WINDOW_SIZE_LDR * 2);
	
	// loop for each pixel in input image
	for (int j = 0; j < (*input).getRows(); j++){		
		for (int i = 0; i < (*input).getCols(); i++){
			numerator = 0.0;			
			denomirator = 0.0;
			weight = 0.0;
			
			// loop for each pixel in kernel
			for (int l = ((j - localWindowSize / 2) > 0) ? (j - localWindowSize / 2) : 0; l <= (j + localWindowSize / 2) && l < (*input).getRows(); l++){
				for (int k = ((i - localWindowSize / 2) > 0) ? (i - localWindowSize / 2) : 0; k <= (i + localWindowSize / 2) && k < (*input).getCols(); k++){
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
 * create grayscale image from RGB image based on color weights / LAB image only using L
 * 
 * @param output - output grayscale image
 * @param input - input image
 */
void TMOBae06::CreateGrayscale(pfstmo::Array2D output, TMOImage * input){
	for (int j = 0; j < input->GetHeight(); j++){		
		for (int i = 0; i < input->GetWidth(); i++){
			// RGB mode
			//output(i, j) = RgbToGray(input->GetPixel(i, j)[0], input->GetPixel(i, j)[1], input->GetPixel(i, j)[2]);			
			
			// LAB mode
			output(i, j) = input->GetPixel(i, j)[0] * SCALE_LAB_TO_RGB;
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
void TMOBae06::GetDetailFromBase(pfstmo::Array2D detail, pfstmo::Array2D base, pfstmo::Array2D input){
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
 * @param height - height of picture in pixels
 * @return computed sigmaS parameter
 */
double TMOBae06::ComputeSigmaS(int width, int height){
	if (hdr){
		// hdr
		return ((height < width) ? height : width) / 64.0;
	}else{	
		// ldr
		return ((height < width) ? height : width) / 16.0;
	}	
}

/**
 * this method performs histogram matching
 * 
 * @param comulativeInputHistogram - normalised cumulative histogram of input picutre (base)
 * @param comulativeModelHistogram - cumulative histogram if model (model base) normalised to same value
 * @param input - array of values for histogram matching, this will be overwritten by new values
 */
void TMOBae06::HistogramMatching(int * comulativeInputHistogram, int * comulativeModelHistogram, pfstmo::Array2D * input){
	int value, inputShade; 		
	
	for (int j = 0; j < (*input).getRows(); j++){
		for (int i = 0; i < (*input).getCols(); i++){
			inputShade = (int) ((*input)(i, j) * HISTOGRAM_BOOSTER_F);
			
			if (inputShade < 0) inputShade = 0;
			if (inputShade > (HISTOGRAM_LEVELS - 1)) inputShade = HISTOGRAM_LEVELS - 1;			
			
			if (verbose) std::cerr << "inputShade: " << inputShade << std::endl;
			
			value = comulativeInputHistogram[inputShade];
			(*input)(i, j) = FindClosestShade(value, comulativeModelHistogram) / HISTOGRAM_BOOSTER_F;
			
			// DEBUG
			/* if (verbose) {
				if ((i % 25 == 0) && (j % 25 == 0)){
					std::cerr << "IN HM inputShade: " << inputShade << std::endl;
					std::cerr << "IN HM value: " << value << std::endl;
					std::cerr << "IN HM result: " << input(i, j) << std::endl << std::endl;
				}
			}*/
		}
	}
}

/**
 * fills the rho variable for textureness transfer
 * 
 * @param rho - rho field, output parameter
 * @param texturenessDesired - T' - input parameter
 * @param texturenessBase - T(B) - input parameter
 * @param texturenessDetail T(D) - input parameter
 */
void TMOBae06::FillRho(pfstmo::Array2D rho, pfstmo::Array2D texturenessDesired, pfstmo::Array2D texturenessBase, pfstmo::Array2D texturenessDetail){
	double max = std::numeric_limits<double>::min();
	
	for (int j = 0; j < texturenessDesired.getRows(); j++){
		for (int i = 0; i < texturenessDesired.getCols(); i++){			
			rho(i, j) = std::max(0.0f, (texturenessDesired(i, j) - texturenessBase(i, j)) / texturenessDetail(i, j));
			
			// find maximum
			if (rho(i, j) > max) max = rho(i, j);
		}
	}
	
	// normalise rho
	if (verbose)  std::cerr << "FillRho max: " << max << std::endl;
	//double divider = (double) max / MAX_RHO;	
	//double divider = 10.0;		
	
	for (int j = 0; j < rho.getRows(); j++){
		for (int i = 0; i < rho.getCols(); i++){
			//rho(i, j) /= divider;
			
			// limit max value of rho to avoid artifacts
			if (rho(i, j) > MAX_RHO) rho(i, j) = MAX_RHO;
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
//void TMOBae06::HighPassFilter(pfstmo::Array2D output, pfstmo::Array2D input, double cutoff){		
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
//}

/**
 * high-pass filter for image processing
 * 
 * @param input - input image, this will be overwritten by filtered image
 */
void TMOBae06::HighPassFilterV2(pfstmo::Array2D * input){		
	pfstmo::Array2D inputBackup = pfstmo::Array2D((*input).getCols(), (*input).getRows());
	
	// create backup array
	for (int j = 0; j < (*input).getRows(); j++){
		for (int i = 0; i < (*input).getCols(); i++){
			inputBackup(i, j) = (*input)(i, j);
		}
	}
	
	// 2D HP filter (http://homepages.inf.ed.ac.uk/rbf/BOOKS/PHILLIPS/cips2ed.pdf p.87)
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
 * computes textureness of given array using cross bilateral filter
 * assumes that sigmaS and sigmaR variables are set
 * 
 * @param texureness - output textureness
 * @param input - input array to compute textureness from
 * @param isDetail - indicate whether or not is this textureness of detail, detail has different scale and doesnt use high-pass filter
 */
void TMOBae06::FillTextureness(pfstmo::Array2D * textureness, pfstmo::Array2D * input, bool isDetail){		
	
	pfstmo::Array2D filteredInput = pfstmo::Array2D((*input).getCols(), (*input).getRows());
	//pfstmo::Array2D* filteredInput = new pfstmo::Array2D(input.getCols(),input.getRows());		
	
	// if (isDetail) -> copy input to filtered input (no need to use HP filter)
	for (int j = 0; j < (*input).getRows(); j++){
		for (int i = 0; i < (*input).getCols(); i++){				
			filteredInput(i, j) = (*input)(i, j) * 256.0;
			//filteredInput(i, j) = (*input)(i, j);
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
	
	// call CBF
	CrossBilateralFilter(textureness, input, &filteredInput);		
}

/**
 * transformation function
 * @return exit code
 */
int TMOBae06::Transform(){
	pSrc->Convert(TMO_LAB);
	pDst->Convert(TMO_LAB);
	
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
	std::string modelFilename = modelFileNameParam.GetString();
	if (verbose) std::cerr << "modelFilename: " << modelFilename << std::endl;
	if (modelFilename == ""){
		std::cerr << "model parameter is mandatory!" << std::endl;
		return -1;
	}

	// load model		
	model = new TMOImage(modelFilename.c_str());	
	model->Convert(TMO_LAB);
	
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
	pfstmo::Array2D rho = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	
	// convert input image and model to grayscale
	CreateGrayscale(grayscale, pSrc);
	CreateGrayscale(modelGrayscale, model);	
	
	// bilateralFiltering for input
	sigmaS = ComputeSigmaS(pSrc->GetWidth(), pSrc->GetHeight());		// proportial to image size
	sigmaR = (hdr) ? SIGMA_R_HDR : SIGMA_R_LDR;
	
	BilateralFilter(&base, &grayscale);
	GetDetailFromBase(detail, base, grayscale);				
	
	// bilateralFiltering for model	
	BilateralFilter(&modelBase, &modelGrayscale);	
	
	// get textureness	
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
	if (verbose){
		PrintHistogram(inputHistogram, "input base histogram");	
		PrintHistogram(modelHistogram, "model base histogram");	
		PrintHistogram(comulativeInputHistogram, "comulative input");
		PrintHistogram(comulativeModelHistogram, "comulative model");
		
		/*PrintHistogram(inputHistogramTextureness, "inputHistogramTextureness");	
		PrintHistogram(modelHistogramTextureness, "modelHistogramTextureness");	
		PrintHistogram(comulativeInputHistogramTextureness, "comulativeInputHistogramTextureness");
		PrintHistogram(comulativeModelHistogramTextureness, "comulativeModelHistogramTextureness");*/
	}
	
	// do historam matching from model base to new input base B'
	HistogramMatching(comulativeInputHistogram, comulativeModelHistogram, &base);
	
	// historam matching for textureness
	HistogramMatching(comulativeInputHistogramTextureness, comulativeModelHistogramTextureness, &texturenessDesired);		
	
	// fill rho array
	FillRho(rho, texturenessDesired, texturenessBase, texturenessDetail);			
	
	// DEBUG
	/*InitialiseHistogram(inputHistogram);
	FillHistogram(base, inputHistogram);
	PrintHistogram(inputHistogram, "input base histogram after matching");	*/
	
	// show debug information
	if (verbose){
		std::cerr << "Rho" << std::endl;
		for (int j = 0; j < rho.getRows(); j+=5){
			for (int i = 0; i < rho.getCols(); i+=5){
				//if (rho(i, j) > 100) std::cerr << "A rho is:" << rho(i, j) << std::endl;
				std::cerr << "rho(" << i << ", " << j << ") " << rho(i, j) << std::endl;
			}
		}
	}

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	double r, g, b, shade;
	int inputShade;
	int value;
	int outputShade;
	
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++){
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++){
			r = *pSourceData++;
			g = *pSourceData++;
			b = *pSourceData++;
			
			// show base
			//shade = base(i, j) / SCALE_LAB_TO_RGB;			
			
			// show detail
			//shade = detail(i, j) / HISTOGRAM_LEVELS_F;
			
			// show base + detail * 3
			//shade = (base(i, j) + detail(i, j) * 3) / HISTOGRAM_LEVELS;
			
			// final result
			shade = (base(i, j) + rho(i, j) * detail(i, j)) / SCALE_LAB_TO_RGB;
			
			// fix "overflows"
			if (shade > 100.0){
				shade = 100.0;
			}else if (shade < 0.0){
				shade = 0.0;
			}
			
			*pDestinationData++ = shade;
			*pDestinationData++ = 0.0;				// creating grayscale in LAB
			*pDestinationData++ = 0.0;
		}
	}		
	
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);		
	
	return 0;
}



