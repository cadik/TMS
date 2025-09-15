/* --------------------------------------------------------------------------- *
 * TMOYu21.h														  		   *
 *																			   *
 * Author: Ludmila Krejčová													   *
 * --------------------------------------------------------------------------- */


#include "TMO.h"
#include <vector>
#include <memory>

class TMOYu21 : public TMO
{
public:
	TMOYu21();
	virtual ~TMOYu21();
	virtual int Transform();
protected:

	using tArray3 = std::array<double, 3>;

	/*
	* Computing constants kr, kg, and kb for the given image based on SSIM (Structural Similarity Index)
	*/
	tArray3 computeK(tArray3 &Krg_Kgb_Kbr, cv::Mat imgR, cv::Mat imgG, cv::Mat imgB);

	/**
 	*
 	* This function computes a contrast-enhanced image using correlation coefficients Krg, Kgb, and Kbr.
 	*
 	* The contrast image is computed (according to the provided sourcecode, not article) as:
 	* \f[
 	* C(x, y) = \cdot (K_{rg} \cdot (R + G) + K_{gb} \cdot (G + B) + K_{br} \cdot (B + R))
 	* \f]
 	*/
	std::unique_ptr<std::vector<double>>  createContrastImage(tArray3 &Krg_Kgb_Kbr);

	/*
	*	This function calculates correlation coefficients (Krg, Kgb, Kbr) 
	*  for the red, green, and blue channels of the input image.
	*/
	tArray3 computeKrg_Kgb_Kbr();


	/*
	* Computing SSIM(R,C)(G,C)(B,C) for this picture
	*/
	tArray3 computeSSIM(std::unique_ptr<std::vector<double>> &contrastImage, cv::Mat imgR, cv::Mat imgG, cv::Mat imgB);

	/*
	* Function to retrieve pixel value for the specified channel (R, G, or B) from the image data
	*/
	inline double getPixel(const double* data, int width, int x, int y, int channel);

	/*
	* Function to set the pixel value for a specific channel (R, G, or B) in the image data
	*/
	inline void setPixel(double* data, int width, int x, int y, int channel, double value);

	/*
	* Function to resize an image using bilinear interpolation
	*/
	std::unique_ptr<double[]> resizeImage(double* input, int srcWidth, int srcHeight, int destWidth, int destHeight);

	/*
	* Compuing intenstity difference from random and neighboring pairs of pixels		
	*/
	std::shared_ptr<std::vector<double>> computeContrastDifferences(const std::vector<std::pair<int, int>> &pairs, const double* image64, 
				const double* image32, int channel);
	
	/*
	* Computes the optimal weights for red, green, and blue channels (wr, wg, wb) that minimize the color energy
	*/
	tArray3 computeWeights(const std::vector<double> &allIr, const std::vector<double> &allIg,
				const std::vector<double> &allIb, const tArray3 &kr_kg_kb);

	
	/* 
	* Computes the color energy of an image based on given weights and a scaling factor.
	* The function's result is highly influenced by the epsilon variable,  
	*  despite the authors claiming otherwise in their paper.
	* 
	* - `w` : Array of three weights (for R, G, B channels).
	* - `k` : Computed constant
	* - `I` : Array of three vectors representing pixel contrasts between color (R, G or B) and contrast image.
	* - `colorIndex` : Index (0 = R, 1 = G, 2 = B) of the color channel used for contrast computation.
	* - `epsilon` : key constant for conversion. Differes for datasets - Cadik 0.15, Color250 0.2, CSDD 0.26  
	*/	
	double 	computeColorEnergy(const tArray3 &w, double k, const std::array<std::vector<double>, 3> &I, size_t colorIndex);

	/*
	* Function to generate random pairs of coordinates within the given size
	*/
	std::vector<std::pair<int, int>> findRandomPairs(int size) const;

	/*
	* Finds if range is 0-1 or in 0-255
	*/
	bool isInRange0to1(double *pSourceData, int numPix);

protected:
	TMODouble epsilonParameter;
	TMOBool HDRParameter;
};
