/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                               *
*                                                                              *
*                       Diploma thesis                                         *
*                       Author: Vladimir Vlkovic                               *
*                       Brno 2017                                              *
*                                                                              *
*                       Implementation of the TMOHu14 class                    *
*                                                                              *
*******************************************************************************/
/**
 * @file TMOHu14.cpp
 * @brief Implementation of the TMOHu14 class
 * @author Vladimir Vlkovic
 * @class TMOHu14
 */

#include "TMOHu14.h"
#include <boost/concept_check.hpp>
#include <complex>

/**
 * @brief Constructor
 */
TMOHu14::TMOHu14()
{
  SetName(L"Hu14");
  SetDescription(L"Image and video decolorization in real-time using dominant colors");
}

/**
 * @brief Destructor
 */
TMOHu14::~TMOHu14()
{
}

/**
 * @brief returns edge mat using Canny edge detector
 * 
 * @param color channel mat
 * @return edge mat of input channel
 */
cv::Mat TMOHu14::getEdgeMat(cv::Mat channel)
{
  double min, max;
  cv::Mat tmp, result;

  //cv::Size size(channel.rows/2,channel.cols/2);
  //cv::Mat data2;
  // cv::resize(src,data2,size);
  int threshold1 = 70;
  int threshold2 = 50;
  cv::minMaxLoc(channel, &min, &max);
  if (min != max)
  {
    channel.convertTo(tmp, CV_8U, 255 * (max - min));
  }
  cv::Canny(tmp, result, threshold1, threshold2);
  if (result.empty())
    result = cv::Mat::zeros(channel.rows, channel.cols, CV_8U); /** in case of frame has no edges e.g. black frame */

  return result;
}

/**
 * @brief Quantizes the colors of the I_edge 
 * 
 * @param src I_edge Mat
 * @param dst quantized I_edge Mat
 */
void TMOHu14::kmeansColorQuantization(const cv::Mat3b &src, cv::Mat3b &dst)
{
  int K = 256; ///should be thi but takes too long

  int n = src.rows * src.cols;
  cv::Mat data = src.reshape(1, n);
  data.convertTo(data, CV_32F);

  std::vector<int> labels;
  cv::Mat1f colors;
  cv::kmeans(data, K, labels, cv::TermCriteria(), 1, cv::KMEANS_PP_CENTERS, colors);

  for (int i = 0; i < n; ++i)
  {
    data.at<float>(i, 0) = colors(labels[i], 0);
    data.at<float>(i, 1) = colors(labels[i], 1);
    data.at<float>(i, 2) = colors(labels[i], 2);
  }

  cv::Mat reduced = data.reshape(3, src.rows);
  reduced.convertTo(dst, CV_8U);

  //cv::imshow("Reduced", dst);
  //  cv::waitKey();
}

/**
 * @brief Converts XYZ vector to bgr
 * http://www.easyrgb.com/en/math.php
 * 
 * @param xyzVector XYZ color vector
 * @return BGR color vector 
 */
cv::Vec3d TMOHu14::xyz2bgr(cv::Vec3d xyzVector)
{
  cv::Vec3d bgrVector;

  double var_X = xyzVector[0] / 100;
  double var_Y = xyzVector[1] / 100;
  double var_Z = xyzVector[2] / 100;

  double r = var_X * 3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
  double g = var_X * -0.9689 + var_Y * 1.8758 + var_Z * 0.0415;
  double b = var_X * 0.0557 + var_Y * -0.2040 + var_Z * 1.0570;

  if (r > 0.0031308)
    r = 1.055 * (pow(r, (1 / 2.4))) - 0.055;
  else
    r = 12.92 * r;
  if (g > 0.0031308)
    g = 1.055 * (pow(g, (1 / 2.4))) - 0.055;
  else
    g = 12.92 * g;
  if (b > 0.0031308)
    b = 1.055 * (pow(b, (1 / 2.4))) - 0.055;
  else
    b = 12.92 * b;

  bgrVector[2] = r * 255;
  bgrVector[1] = g * 255;
  bgrVector[0] = b * 255;

  return bgrVector;
}

/**
 * @brief Converts Luv vector to bgr
 * http://www.easyrgb.com/en/math.php
 * 
 * @param luvVector Luv color vector
 * @return BGR color vector 
 */
cv::Vec3d TMOHu14::Luv2rgb(cv::Vec3d luvVector)
{
  cv::Vec3d bgrVector, xyzVector;
  double Y_n = 100.0;
  double u_line_n = 0.2009;
  double v_line_n = 0.4610;

  double L = luvVector[0];
  double u = luvVector[1];
  double v = luvVector[2];

  double ref_X = 95.047;
  double ref_Y = 100.0;
  double ref_Z = 108.883;
  double x, y, z;

  double var_Y = (L + 16) / 116;
  if (pow(var_Y, 3) > 0.008856)
    var_Y = pow(var_Y, 3);
  else
    var_Y = (var_Y - 16 / 116) / 7.787;

  double ref_U = (4 * ref_X) / (ref_X + (15 * ref_Y) + (3 * ref_Z));
  double ref_V = (9 * ref_Y) / (ref_X + (15 * ref_Y) + (3 * ref_Z));

  double var_U = u / (13 * L) + ref_U;
  double var_V = v / (13 * L) + ref_V;

  y = var_Y * 100;
  x = -(9 * y * var_U) / ((var_U - 4) * var_V - var_U * var_V);
  z = (9 * y - (15 * var_V * y) - (var_V * x)) / (3 * var_V);

  xyzVector[0] = x;
  xyzVector[1] = y;
  xyzVector[2] = z;

  bgrVector = TMOHu14::xyz2bgr(xyzVector);

  return bgrVector;
}

/**
 * @brief Converts BGR vector to XYZ and then to Luv
 * @param bgrVector BGR color vector
 * @return Luv color vector 
 */
cv::Vec3d TMOHu14::rgb2Luv(cv::Vec3i bgrVector)
{
  double x, y, z, L, u, v, r, g, b;
  cv::Vec3d LuvVector;

  r = static_cast<double>(bgrVector[2]) / 255.0f;
  g = static_cast<double>(bgrVector[1]) / 255.0f;
  b = static_cast<double>(bgrVector[0]) / 255.0f;

  if (r > 0.04045)
    r = std::pow((r + 0.055) / 1.055, 2.4);
  else
    r = r / 12.92;

  if (g > 0.04045)
    g = std::pow((g + 0.055) / 1.055, 2.4);
  else
    g = g / 12.92;

  if (b > 0.04045)
    b = std::pow((b + 0.055) / 1.055, 2.4);
  else
    b = b / 12.92;

  r = r * 100;
  g = g * 100;
  b = b * 100;

  x = r * 0.4124 + g * 0.3576 + b * 0.1805;
  y = r * 0.2126 + g * 0.7152 + b * 0.0722;
  z = r * 0.0193 + g * 0.1192 + b * 0.9505;

  TMOImage::XyzToLuv(x, y, z, &L, &u, &v);

  LuvVector[0] = L;
  LuvVector[1] = u;
  LuvVector[2] = v;

  return LuvVector;
}

/**
 * @brief Get feature vector consisting of Luv color and color percentage in image
 * @param src I_edge Mat
 * @return feature vector/color pallete/color histogram : color in Luv, color percentage in image
 */
void TMOHu14::getPalette(std::map<cv::Vec3d, float, lessVec3b> &paletteLuv, cv::Mat &src)
{
  std::map<cv::Vec3i, float, lessVec3b> paletteRGB;
  // std::map<cv::Vec3d, int, lessVec3b> paletteLuv;
  float pixelCount = 0.0; //src.rows*src.cols;
  for (int r = 0; r < src.rows; ++r)
  {
    for (int c = 0; c < src.cols; ++c) /** get every color and pixel count of every color*/
    {

      cv::Vec3b color = src.at<cv::Vec3b>(r, c);

      if (!(color[0] == 0 && color[1] == 0 && color[2] == 0))
      {
        if (paletteRGB.count(color) == 0)
        {
          paletteRGB[color] = 1;
        }
        else
        {
          paletteRGB[color] = paletteRGB[color] + 1;
        }
        pixelCount++;
      }
    }
  }
  int d1e = paletteRGB.size();
  std::map<cv::Vec3i, float, lessVec3b>::iterator it = paletteRGB.begin();
  while (it != paletteRGB.end())
  {
    float colorPercentage = 0.0;
    colorPercentage = (it->second / pixelCount);
    if (colorPercentage < 0.001f) /** if color percentage is less then discard */
    {
      paletteRGB.erase(it++);
    }
    else
    {
      cv::Vec3d tmpBgr = TMOHu14::rgb2Luv(it->first); /**convert bgr to Luv*/
      paletteLuv[tmpBgr] = colorPercentage;
      ++it;
    }
  }
}

/**
 * @brief Get get dominant color descriptr in LUV space merging of perceptual similar colors
 * 
 * @param palette Luv color pallete with pixel percentage
 * @return dominat color feature vector : color in Luv, color percentage in image
 */
std::map<cv::Vec3d, float, lessVec3b> TMOHu14::getDominantColorDescriptor(std::map<cv::Vec3d, float, lessVec3b> palette)
{

  int d = 0;
  cv::Vec3d newColor;
  std::map<cv::Vec3d, float, lessVec3b> palette2;
  std::map<cv::Vec3d, float, lessVec3b>::iterator it = palette.begin();
  std::map<cv::Vec3d, float, lessVec3b>::iterator it2 = palette.begin();
  it++;
  while (it != palette.end())
  {
    double L, u, v;
    it2 = palette.begin();

    L = it->first[0];
    u = it->first[1];
    v = it->first[2];

    while (it2 != palette.end())
    {
      double L2, u2, v2;
      L2 = it2->first[0];
      u2 = it2->first[1];
      v2 = it2->first[2];
      d = std::sqrt(std::pow(it->first[0] - L2, 2) + std::pow(it->first[1] - u2, 2) + std::pow(it->first[2] - v2, 2));

      if (d < DELTA && d > 0)
      {

        newColor[0] = (it->first[0] * it->second + L2 * it2->second) / (it->second + it2->second);
        newColor[1] = (it->first[1] * it->second + u2 * it2->second) / (it->second + it2->second);
        newColor[2] = (it->first[2] * it->second + v2 * it2->second) / (it->second + it2->second);

        palette.erase(it++);
        palette.erase(it2++);
        palette[newColor] = it->second + it2->second;

        break;
      }
      else
        ++it2;
    }
    ++it;
  }

  return palette;
}

/**
 * @brief Get grayscale pallete form bgr pallete using weights 
 * 
 * @param weight_r weight for red color in range 0.0-1.0
 * @param weight_g weight for green color in range 0.0-1.0
 * @param weight_b weight for blue color in range 0.0-1.0
 * @param luvBgrPalette vector LUV and color percentage[3], second BGR
 * @return first: luv color and percentage, second: grayscale value
 */
std::map<cv::Vec4d, int, lessVec4d> TMOHu14::getGrayscalePalette(float weight_r, float weight_g, float weight_b, std::map<cv::Vec4d, cv::Vec3d, lessVec4d> luvBgrPalette)
{
  int grayscaleColor = 0;
  std::map<cv::Vec4d, int, lessVec4d> grayscalePalette;

  for (std::map<cv::Vec4d, cv::Vec3d, lessVec4d>::iterator it = luvBgrPalette.begin(); it != luvBgrPalette.end(); ++it)
  {
    grayscaleColor = it->second[2] * weight_r + it->second[1] * weight_g + it->second[0] * weight_b;
    grayscalePalette[it->first] = grayscaleColor;
  }

  return grayscalePalette;
}
/**
 * @brief Get the H metric 
 * 
 * @param color1 first colo ro be used in the computation
 * @param color2 second color to be used in the computation
 * @return the computed H metric
 **/
double TMOHu14::getHMetric(cv::Vec4d color1, cv::Vec4d color2)
{

  double dist, h;

  dist = std::sqrt(std::pow(color1[0] - color2[0], 2) + std::pow(color1[1] - color2[1], 2) + std::pow(color1[2] - color2[2], 2));
  h = dist * (color1[3] + color2[3]);
  return h;
}
/**
 * @brief Get the Xi metric which determines which combination of weights to use
 * @param grayscalePalette  first: luv color and percentage, second: grayscale value
 * @return Xi metric
 */
double TMOHu14::getXiMetric(std::map<cv::Vec4d, int, lessVec4d> grayscalePalette)
{

  double h = 0.0;
  float k = 0.0;
  float d = 0.0;
  float dist = 0.0;
  double xi = 0.0;
  int i = 0;

  for (std::map<cv::Vec4d, int, lessVec4d>::iterator itGray1 = grayscalePalette.begin(); itGray1 != grayscalePalette.end(); ++itGray1)
  {
    for (std::map<cv::Vec4d, int, lessVec4d>::iterator itGray2 = grayscalePalette.begin(); itGray2 != grayscalePalette.end(); ++itGray2)
    {
      h = TMOHu14::getHMetric(itGray1->first, itGray2->first);
      d = abs(itGray1->second - itGray2->second);
      if (d > TAU)
        k = 1;
      else
        k = 0;
      xi += h * (LAMBDA * k + (1 - LAMBDA) * d);
      i++;
    }
  }
  return xi;
}
/**
 * @brief Sorting function for sorting the 4d vector in descending order by the first vector member
 */
bool TMOHu14::sortFunc(const cv::Vec4d &a, const cv::Vec4d &b) { return a[0] > b[0]; }

/**
 * @brief Get the best weights combination according to the Xi metric
 * 
 * @param luvBgrPalette vector LUV and color percentage[3], second BGR
 * @return vector of ordered weights from best to worst
 */
std::vector<cv::Vec4d> TMOHu14::getBestWeightsCandidate(std::map<cv::Vec4d, cv::Vec3d, lessVec4d> luvBgrPalette, cv::Mat redMat, cv::Mat greenMat, cv::Mat blueMat)
{
  std::vector<cv::Vec4d> weights;
  cv::Vec4d data;
  float weight_r = 0;
  float weight_g = 0;
  float weight_b = 0;
  double maxXi = 0;
  double xi;

  std::map<cv::Vec4d, int, lessVec4d> grayscalePalette;

  for (weight_r = 0; weight_r <= 10; weight_r++)
  {
    for (weight_g = 0; weight_g <= 10; weight_g++)
    {
      weight_b = 10 - (weight_r + weight_g);
      if (weight_b >= 0)
      {
        float weight_r_f = weight_r / 10.0;
        float weight_g_f = weight_g / 10.0;
        float weight_b_f = weight_b / 10.0;
        grayscalePalette = TMOHu14::getGrayscalePalette(weight_r_f, weight_g_f, weight_b_f, luvBgrPalette);
        xi = TMOHu14::getXiMetric(grayscalePalette);

        //cv::imshow("Reduced", redMat * weight_r_f + greenMat * weight_g_f + blueMat*weight_b_f);
        //cv::waitKey();
        //if(xi > maxXi)
        //{
        data[0] = xi;

        data[1] = weight_r_f;
        data[2] = weight_g_f;
        data[3] = weight_b_f;
        //}
        //count++;
        weights.push_back(data);
      }
    }
  }
  std::sort(weights.begin(), weights.end(), sortFunc);

  return weights;
}

/**
 * @brief Gets the histogram for the frame
 * 
 * @param frame video frame which histogram we want
 * @return histogram of the frame
 * */
cv::Mat TMOHu14::getHistogram(cv::Mat frame)
{
  cv::Mat luvFrame;
  int lbins = 30, ubins = 30, vbins = 30;
  float lrange[] = {0, 256};
  float urange[] = {0, 256};
  float vrange[] = {0, 256};
  const float *ranges[] = {lrange, urange, vrange};
  int histSize[] = {lbins, ubins, vbins};
  cv::Mat hist;
  int channelss[] = {0, 1, 2};
  cv::cvtColor(frame, luvFrame, cv::COLOR_BGR2Luv, 0);
  cv::normalize(luvFrame, luvFrame, 0.0, 255.0, cv::NORM_MINMAX, CV_32F);

  cv::calcHist(&luvFrame, 1, channelss, cv::Mat(), hist, 3, histSize, ranges, true, false);
  cv::normalize(hist, hist, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());

  return hist;
}
/**
 * @brief Gets the changig rate between two histograms
 * 
 * @param histogramDifference the difference between currect frame histogram and previous frame histogram
 * @return changing rate
 * */
double TMOHu14::getChangingRate(double histogramDifference)
{
  if (histogramDifference < THETA)
  {
    return std::pow(THETA * (histogramDifference / THETA), GAMMA);
  }
  else
  {
    double tmp;
    tmp = (1 - histogramDifference) / (1 - THETA);
    return 1 + (THETA - 1) * std::pow(tmp, GAMMA);
  }
}
/**
 * @brief Returns the closest weight to the weight of the previous frame
 * 
 * @param previousWeight the weight that was used to decolorize the previous frame
 * @param  currentWeights list of ordered weights that can be considered to decolorize current frame
 * @param searchNumber indicates how deep in the ordered weights should be searched
 * @return vector of the closest weight
 * */
cv::Vec4d TMOHu14::getClosestWeight(cv::Vec4d previousWeight, std::vector<cv::Vec4d> currentWeights, int searchNumber)
{
  double prevRedWeight, prevGreenWeight, prevBlueWeight, currRedWeight, currGreenWeight, currBlueWeight;
  double proximity;
  double minProx;
  int res = 0;
  prevRedWeight = previousWeight[1];
  prevGreenWeight = previousWeight[2];
  prevBlueWeight = previousWeight[3];
  for (int i = 0; i < searchNumber; i++)
  {
    currRedWeight = currentWeights[i][1];
    currGreenWeight = currentWeights[i][2];
    currBlueWeight = currentWeights[i][3];

    proximity = std::sqrt(std::pow(prevRedWeight - currRedWeight, 2) + std::pow(prevGreenWeight - currGreenWeight, 2) + std::pow(prevBlueWeight - currBlueWeight, 2));
    if (i == 0)
      minProx = proximity;
    if (proximity < minProx)
      res = i;
  }
  return currentWeights[res];
}

/**
 * @brief High-Performance Flexible Broadband Photodetector Based on Organolead Halide Perovskite
 */
int TMOHu14::Transform()
{
  double *pSourceData = pSrc->GetData();
  double *pDestinationData = pDst->GetData();

  int height = pSrc->GetHeight();
  int width = pSrc->GetWidth();

  double min, max;
  cv::Mat redMat, greenMat, blueMat, redEdgeMat, greenEdgeMat, blueEdgeMat, sumEdgeMat, tmpMat, r_down, g_down, b_down;
  /** Mat for each color channel */
  cv::Mat3b reduced;
  double histogramDifference = 0.0, changingRate = 0.0;
  int heightDown = height / 2;
  int widthDown = width / 2;
  cv::Mat mergedMat;
  cv::Size size(width / 2, height / 2);

  redMat = cv::Mat::zeros(height, width, CV_32F);
  greenMat = cv::Mat::zeros(height, width, CV_32F); /** mats for color channels*/
  blueMat = cv::Mat::zeros(height, width, CV_32F);

  redEdgeMat = cv::Mat::zeros(heightDown, widthDown, CV_8U);
  greenEdgeMat = cv::Mat::zeros(heightDown, widthDown, CV_8U); /** mats for edges of color channels must be CV_8U because Canny function */
  blueEdgeMat = cv::Mat::zeros(heightDown, widthDown, CV_8U);

  std::map<cv::Vec3d, float, lessVec3b> palette;
  std::map<cv::Vec4d, cv::Vec3d, lessVec4d> luvBgrPalette; /**first vector LUV and color percentage, second BGR*/

  std::vector<cv::Vec4d> weights; /** vector of 66 weights, weight format is [0]the xi number [1]red weight [2] green weight [3]blue weight */

  for (int j = 0; j < pSrc->GetHeight(); j++)
  {
    pSrc->ProgressBar(j, pSrc->GetHeight());
    for (int i = 0; i < pSrc->GetWidth(); i++)
    {
      redMat.at<float>(j, i) = *pSourceData++;
      greenMat.at<float>(j, i) = *pSourceData++; /** getting separate RGB channels */
      blueMat.at<float>(j, i) = *pSourceData++;
    }
  }

  cv::normalize(redMat, redMat, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);
  cv::normalize(greenMat, greenMat, 0.0, 1.0, cv::NORM_MINMAX, CV_32F); /**normalize */
  cv::normalize(blueMat, blueMat, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);

  cv::resize(redMat, r_down, size);
  cv::resize(greenMat, g_down, size); /** img downsamplig */
  cv::resize(blueMat, b_down, size);

  redEdgeMat = TMOHu14::getEdgeMat(r_down);
  greenEdgeMat = TMOHu14::getEdgeMat(g_down); /** getting edge mats */
  blueEdgeMat = TMOHu14::getEdgeMat(b_down);

  sumEdgeMat = redEdgeMat + greenEdgeMat + blueEdgeMat; /** need to sum all channel edges */
  sumEdgeMat.convertTo(tmpMat, CV_32F);                 /** conversion to float mat*/

  redEdgeMat = r_down.mul(tmpMat);
  greenEdgeMat = g_down.mul(tmpMat); /** multipling edge map by color channel maps to achieve I_edge (see alg.) */
  blueEdgeMat = b_down.mul(tmpMat);

  std::vector<cv::Mat> channels;
  channels.push_back(blueEdgeMat);
  channels.push_back(greenEdgeMat);
  channels.push_back(redEdgeMat); /** mergig channels into one mat */
  cv::merge(channels, mergedMat);
  cv::normalize(mergedMat, mergedMat, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);

  if (cv::countNonZero(blueEdgeMat) > 0 && cv::countNonZero(greenEdgeMat) > 0 && cv::countNonZero(redEdgeMat) > 0) /** if there are edges in image do */
  {
    TMOHu14::kmeansColorQuantization(mergedMat, reduced); /** decrease number of colors uing color quantization */

    TMOHu14::getPalette(palette, reduced); /** gettign the color palette in luv */

    palette = TMOHu14::getDominantColorDescriptor(palette);
    int e = palette.size();
    while (1)
    {
      palette = TMOHu14::getDominantColorDescriptor(palette); /** cycle util palette cannot be smaller */
      if (e == palette.size())
      {
        break;
      }
      e = palette.size();
    }

    for (std::map<cv::Vec3d, float, lessVec3b>::iterator it = palette.begin(); it != palette.end(); ++it)
    {
      cv::Vec3d tmp;
      cv::Vec4d tmp2;
      tmp = TMOHu14::Luv2rgb(it->first);
      tmp2[0] = it->first[0];
      tmp2[1] = it->first[1];
      tmp2[2] = it->first[2];
      tmp2[3] = it->second;
      luvBgrPalette[tmp2] = tmp; /** creation of first vector LUV and color percentage, second BGR */
    }

    weights = getBestWeightsCandidate(luvBgrPalette, redMat, greenMat, blueMat); /** gets the ordered vector of weights from best to worst */
  }
  else
  {
    cv::Vec4d tempWeight = 0;
    weights.push_back(tempWeight);
  }

  for (int j = 0; j < pSrc->GetHeight(); j++)
  {
    pSrc->ProgressBar(j, pSrc->GetHeight());

    for (int i = 0; i < pSrc->GetWidth(); i++) /** result to output, taking only the image correction is discarded */
    {
      double final = redMat.at<float>(j, i) * weights[0][1] + greenMat.at<float>(j, i) * weights[0][2] + blueMat.at<float>(j, i) * weights[0][3];

      *pDestinationData++ = final;
      *pDestinationData++ = final;
      *pDestinationData++ = final;
    }
  }

  return 0;
}

int TMOHu14::TransformVideo()
{
  int height = vSrc->GetHeight();
  int width = vSrc->GetWidth();
  cv::VideoCapture cap = vSrc->getVideoCaptureObject();

  double min, max;
  cv::Mat redMat, greenMat, blueMat, redEdgeMat, greenEdgeMat, blueEdgeMat, sumEdgeMat, tmpMat, r_down, g_down, b_down, frame, frame2, previousHist,
      currentHist; /** Mat for each color channel */
  cv::Mat3b reduced;
  cv::Mat bgr[3];
  double histogramDifference = 0.0, changingRate = 0.0;
  int heightDown = height / 2;
  int widthDown = width / 2;
  int searchNumber = 22;
  cv::Mat mergedMat;
  cv::Size size(width / 2, height / 2);

  redMat = cv::Mat::zeros(height, width, CV_32F);
  greenMat = cv::Mat::zeros(height, width, CV_32F); /** mats for color channels */
  blueMat = cv::Mat::zeros(height, width, CV_32F);

  redEdgeMat = cv::Mat::zeros(heightDown, widthDown, CV_8U);
  greenEdgeMat = cv::Mat::zeros(heightDown, widthDown, CV_8U); /** mats for edges of color channels must be CV_8U because Canny function */
  blueEdgeMat = cv::Mat::zeros(heightDown, widthDown, CV_8U);

  std::map<cv::Vec3d, float, lessVec3b> palette;
  std::map<cv::Vec4d, cv::Vec3d, lessVec4d> luvBgrPalette; /** first vector LUV and color percentage, second BGR */

  cv::Mat temp, result;
  std::vector<cv::Vec4d> currentWeights; /** vector of 66 weights, weight format is [0]the xi number [1]red weight [2] green weight [3]blue weight */
  cv::Vec4d previousWeight, closestWeight, finalCurrWeight;
  ; /** look above */

  for (int frameCounter = 0; frameCounter < vSrc->GetTotalNumberOfFrames(); frameCounter++)
  {
    vSrc->GetMatVideoFrame(cap, frameCounter, frame); /** read frame */

    cv::split(frame, bgr);

    cv::normalize(bgr[2], redMat, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);
    cv::normalize(bgr[1], greenMat, 0.0, 1.0, cv::NORM_MINMAX, CV_32F); /** split into rgb channels */
    cv::normalize(bgr[0], blueMat, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);

    cv::resize(redMat, r_down, size);
    cv::resize(greenMat, g_down, size); /** img downsamplig */
    cv::resize(blueMat, b_down, size);

    redEdgeMat = TMOHu14::getEdgeMat(r_down);
    greenEdgeMat = TMOHu14::getEdgeMat(g_down); /** getting edge mats */
    blueEdgeMat = TMOHu14::getEdgeMat(b_down);

    sumEdgeMat = redEdgeMat + greenEdgeMat + blueEdgeMat; /** need to sum all channel edges */
    sumEdgeMat.convertTo(tmpMat, CV_32F);                 /** conversion to float mat */

    redEdgeMat = r_down.mul(tmpMat);
    greenEdgeMat = g_down.mul(tmpMat); /** multipling edge map by color channel maps to achieve I_edge (see alg.) */
    blueEdgeMat = b_down.mul(tmpMat);

    std::vector<cv::Mat> channels;
    channels.push_back(blueEdgeMat);
    channels.push_back(greenEdgeMat);
    channels.push_back(redEdgeMat); /** mergig channels into one mat */
    cv::merge(channels, mergedMat);
    cv::normalize(mergedMat, mergedMat, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);

    if (cv::countNonZero(blueEdgeMat) > 0 && cv::countNonZero(greenEdgeMat) > 0 && cv::countNonZero(redEdgeMat) > 0) /** if there are edges in image do */
    {
      TMOHu14::kmeansColorQuantization(mergedMat, reduced); /** decrease number of colors uing color quantization */

      TMOHu14::getPalette(palette, reduced); /** gettign the color palette in luv */

      palette = TMOHu14::getDominantColorDescriptor(palette);
      int e = palette.size();
      while (1)
      {
        palette = TMOHu14::getDominantColorDescriptor(palette); /** cycle util palette cannot be smaller */
        if (e == palette.size())
        {
          break;
        }
        e = palette.size();
      }

      for (std::map<cv::Vec3d, float, lessVec3b>::iterator it = palette.begin(); it != palette.end(); ++it)
      {
        cv::Vec3d tmp;
        cv::Vec4d tmp2;
        tmp = TMOHu14::Luv2rgb(it->first);
        tmp2[0] = it->first[0];
        tmp2[1] = it->first[1];
        tmp2[2] = it->first[2];
        tmp2[3] = it->second;
        luvBgrPalette[tmp2] = tmp; /** creation of first vector LUV and color percentage, second BGR */
      }

      currentWeights = getBestWeightsCandidate(luvBgrPalette, redMat, greenMat, blueMat); /** gets the ordered vector of weights from best to worst */
    }
    else
    {
      cv::Vec4d tempWeight = 0;
      currentWeights.push_back(tempWeight);
    }

    currentHist = getHistogram(frame); /** get the histogram of the current frame */

    if (frameCounter != 0) /** if its not the first frame */
    {
      histogramDifference = cv::compareHist(currentHist, previousHist, cv::HISTCMP_BHATTACHARYYA); /** compute histogram diff */
      changingRate = getChangingRate(histogramDifference);                                     /** changig rate */

      if (histogramDifference > THETA)
        searchNumber = 2; /** set search number */

      closestWeight = getClosestWeight(previousWeight, currentWeights, searchNumber); /** closest weight */
      finalCurrWeight[1] = previousWeight[1] + changingRate * (closestWeight[1] - previousWeight[1]);
      finalCurrWeight[2] = previousWeight[2] + changingRate * (closestWeight[2] - previousWeight[2]); /** compute final weights */
      finalCurrWeight[3] = previousWeight[3] + changingRate * (closestWeight[3] - previousWeight[3]);
      previousWeight = finalCurrWeight; /**set this current weight to previous for another iteration */
    }
    else /** if frame is first */
    {

      finalCurrWeight = currentWeights[0]; /** set the best weight as final */
      previousWeight = currentWeights[0];  /** set the best weight as previous for another iteration */
      previousHist = currentHist;          /** set histogram to previous for another iteration */
    }

    temp = cv::Mat::zeros(height, width, CV_64FC1);
    result = cv::Mat::zeros(height, width, CV_64FC3);
    for (int j = 0; j < vSrc->GetHeight(); j++)
    {
      for (int i = 0; i < vSrc->GetWidth(); i++) /** result to output, taking only the image correction is discarded */
      {

        temp.at<double>(j, i) = redMat.at<float>(j, i) * finalCurrWeight[1] + greenMat.at<float>(j, i) * finalCurrWeight[2] + blueMat.at<float>(j, i) * finalCurrWeight[3]; /** set the colors */
      }
    }
    channels.clear();
    channels.push_back(temp);
    channels.push_back(temp);
    channels.push_back(temp);
    cv::merge(channels, result); /** create matrix */

    std::cerr << frameCounter + 1 << " / " << vSrc->GetTotalNumberOfFrames() << std::endl;
    vDst->setMatFrame(vDst->getVideoWriterObject(), result); /** write frame; */
  }
  return 0;
}
