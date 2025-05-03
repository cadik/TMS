/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                 A tone reproduction operator for all luminance               *
 *                   ranges considering human color perception                  *
 * 																			                                        *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2024                                 *
 *                                                                              *
 *******************************************************************************/

#include "TMOMikamo14.h"

TMOMikamo14::TMOMikamo14()
{
  SetName(L"Mikamo14");
  SetDescription(L"A tone reproduction operator for all luminance ranges considering human color perception. Two optional parameters, if both set, just ari is used.");

  ari.SetName(L"ari");
  ari.SetDescription(L"Adapted retinal illuminance (ari) in Trolands; <0.0, 1000.0> (optional)");
  ari.SetDefault(10.0);
  ari = 10.0;
  ari.SetRange(0.0, 1000.0);
  this->Register(ari);

  al.SetName(L"al");
  al.SetDescription(L"Adapted luminance (al) in cd/m^2; <0.0, 1000.0> (optional)");
  al.SetDefault(0.0);
  al = 0.0;
  al.SetRange(0.0, 1000.0);
  this->Register(al);

  nob.SetName(L"nob");
  nob.SetDescription(L"Number of bins to be used in wavelength discrimination. Higher number, more accurate. (optional)");
  nob.SetDefault(120);
  nob = 120;
  nob.SetRange(10, 200);
  this->Register(nob);

  rcf.SetName(L"rcf");
  rcf.SetDescription(L"Range compression factor. Amount of compression should be based on type of vision. <0.25, 1.0> (optional)");
  rcf.SetDefault(0.5);
  rcf = 0.5;
  rcf.SetRange(0.25, 1.0);
  this->Register(rcf);
}

TMOMikamo14::~TMOMikamo14() {}

double **TMOMikamo14::getNewColorData()
{
  double **newColorData = new double *[nob];

  // for each bin create an array
  for (int i = 0; i < nob; i++)
  {
    newColorData[i] = new double[colors];
    // each array filled by computed data
    for (int j = 0; j < colors; j++)
    {
      // fraction representing position of new bin in the original data
      double fract = double(i) / double(nob - 1) * double(bins - 1);
      int index = std::floor(fract);
      // compute difference to be added to the origanl value compensating the fractional difference in the position of the new bin
      double diff = color_data[std::min(index + 1, bins - 1)][j] - color_data[index][j];

      newColorData[i][j] = color_data[index][j] + diff * (fract - double(index));
    }
  }

  return newColorData;
}

double **TMOMikamo14::getNewLMSSens()
{
  double **newLMSSens = new double *[nob];

  // for each bin create an array
  for (int i = 0; i < nob; i++)
  {
    newLMSSens[i] = new double[3];
    // each array filled by computed data
    for (int j = 0; j < 3; j++)
    {
      // fraction representing position of new bin in the original data
      double fract = double(i) / double(nob - 1) * double(bins - 1);
      int index = std::floor(fract);
      // compute difference to be added to the origanl value compensating the fractional difference in the position of the new bin
      double diff = LMSsensitivities[std::min(index + 1, bins - 1)][j] - LMSsensitivities[index][j];

      newLMSSens[i][j] = LMSsensitivities[index][j] + diff * (fract - double(index));
    }
  }

  return newLMSSens;
}

double TMOMikamo14::getAdaptedRetinalIlluminance()
{
  // if adapted luminance is set, return it multiplied by the pupil area
  if (al != 0.0)
  {
    double pupilDiameter = 5.697 - 0.658 * std::log10(al) + 0.07 * std::pow(std::log10(al), 2); // pupil diameter depending on the adapted luminance, equation by Blackie and Howland (1999)
    double area = M_PI * std::pow(pupilDiameter / 2, 2);                                        // area of the pupil
    return al * area;
  }

  // if adapted luminance is not set, return the adapted retinal illuminance
  return ari;
}

std::vector<double> TMOMikamo14::getDiscriminationParams(double I)
{
  std::vector<double> params(9);

  params[0] = -18.3 / (1 + 7.2 * std::pow(I, -0.7)) - 0.9;   // λl(I)
  params[1] = -44.6 / (1 + 35.4 * std::pow(I, -1.2)) + 22.0; // λm(I)
  params[2] = 43.0 / (1 + 9.0 * std::pow(I, -1.5)) + 28.0;   // λs(I)

  params[3] = 6.69 / (1 + 2500 * std::pow(I, -2.65)) + 0.80; // k1(I)
  params[4] = -6.24 / (1 + 2500 * std::pow(I, -2.5)) - 0.77; // k2(I)
  params[5] = 0.36 / (1 + 50.02 * std::pow(I, -1.5)) + 0.04; // k3(I)

  params[6] = 0.24 / (1 + 50.04 * std::pow(I, -1.7)) + 0.03; // k4(I)
  params[7] = 0.42 / (1 + 1.76 * std::pow(I, -0.02)) + 0.14; // k5(I)
  params[8] = 0.15 / (1 + 2.80 * std::pow(I, -0.46)) - 0.27; // k6(I)

  return params;
}

double TMOMikamo14::lambdaAdjustment(double ***newLMSSens, int i, int step, int cone)
{
  if ((i - step < 0) || (i - step >= nob))
  {
    return 0.0;
  }
  else
  {
    return (*newLMSSens)[i - step][cone];
  }
}

cv::Mat TMOMikamo14::applyTwoStageModel(double ***newLMSSens, std::vector<double> spd, double I, std::vector<double> params)
{
  // initialize opponent color values
  double V = 0.0;
  double Org = 0.0;
  double Oyb = 0.0;

  double newBinWidth = double(upperBound - lowerBound) / double(nob);

  // go through the bins to get integrated opponent color values
  for (int i = 0; i < nob; i++)
  {
    // matrix with horizontally moved spectral sensitivities
    cv::Mat CmClCs = (cv::Mat_<double>(3, 1) << lambdaAdjustment(newLMSSens, i, (int)(params[0] / binWidth), 0), lambdaAdjustment(newLMSSens, i, (int)(params[1] / binWidth), 1), lambdaAdjustment(newLMSSens, i, (int)(params[2] / binWidth), 2));
    // matrix which adjusts the amplitudes of the cone responses
    cv::Mat M = (cv::Mat_<double>(3, 3) << 0.6, 0.4, 0.0, params[3], params[4], params[5], params[6], params[7], params[8]);
    // get spectral opponent color values
    cv::Mat z = M * CmClCs;
    // add the spectral opponent color values to the integrated opponent color values
    V += spd[i] * z.at<double>(0, 0) * newBinWidth;
    Org += spd[i] * z.at<double>(1, 0) * newBinWidth;
    Oyb += spd[i] * z.at<double>(2, 0) * newBinWidth;
  }
  // create matrix with opponent color values
  cv::Mat opponentColor = (cv::Mat_<double>(3, 1) << V, Org, Oyb);
  // adjust the gap in the viewing conditions
  std::vector<double> newParams = getDiscriminationParams(150.0);
  cv::Mat invM = (cv::Mat_<double>(3, 3) << 0.6, 0.4, 0.0, newParams[3], newParams[4], newParams[5], newParams[6], newParams[7], newParams[8]);
  invM = invM.inv();
  opponentColor = invM * opponentColor;

  return opponentColor;
}

std::vector<double> TMOMikamo14::RGBtoSpectrum(double ***newColorData, double red, double green, double blue)
{
  std::vector<double> spectrum;

  for (int binIndex = 0; binIndex < nob; binIndex++)
  {
    // get spectral power distribution for each color in the bin
    double white_spd = (*newColorData)[binIndex][White];
    double cyan_spd = (*newColorData)[binIndex][Cyan];
    double magenta_spd = (*newColorData)[binIndex][Magenta];
    double yellow_spd = (*newColorData)[binIndex][Yellow];
    double red_spd = (*newColorData)[binIndex][Red];
    double green_spd = (*newColorData)[binIndex][Green];
    double blue_spd = (*newColorData)[binIndex][Blue];
    double spd = 0.0;
    // algorithm to convert RGB to spectral power distribution
    if (red <= green && red <= blue)
    {
      spd += white_spd * red;
      if (green <= blue)
      {
        spd += cyan_spd * (green - red);
        spd += blue_spd * (blue - green);
      }
      else
      {
        spd += cyan_spd * (blue - red);
        spd += green_spd * (green - blue);
      }
    }
    else if (green <= red && green <= blue)
    {
      spd += white_spd * green;
      if (red <= blue)
      {
        spd += magenta_spd * (red - green);
        spd += blue_spd * (blue - red);
      }
      else
      {
        spd += magenta_spd * (blue - green);
        spd += red_spd * (red - blue);
      }
    }
    else
    {
      spd += white_spd * blue;
      if (red <= green)
      {
        spd += yellow_spd * (red - blue);
        spd += green_spd * (green - red);
      }
      else
      {
        spd += yellow_spd * (green - blue);
        spd += red_spd * (red - green);
      }
    }
    // add the spectral power distribution to the spectrum
    spectrum.push_back(spd);
  }
  return spectrum;
}

double TMOMikamo14::luminanceReduction(double Y, double YLogAvg, double Ymax)
{
  // get key value for luminance reduction
  double alpha = 1.03 - 2 / (2 + std::log10(YLogAvg + 1));
  // compute reduced luminance
  double Yr = (alpha * Y) / YLogAvg;
  // compute final, normalized luminance
  double Yn = (Yr * (1 + (Yr / std::pow(Ymax, 2)))) / (1 + Yr);
  return Yn;
}

int TMOMikamo14::Transform()
{
  // adjust number of used bins to change accuracy
  double **newColorData = getNewColorData();
  double **newLMSSens = getNewLMSSens();

  double *pSourceData = pSrc->GetData();
  double *pDestinationData = pDst->GetData();
  double I = getAdaptedRetinalIlluminance();

  // get discrimination parameters for given adapted retinal illuminance
  std::vector<double> params = getDiscriminationParams(I);

  // go through the image and apply the tone mapping operator
  for (int y = 0; y < pSrc->GetHeight(); y++)
  {
    for (int x = 0; x < pSrc->GetWidth(); x++)
    {
      double *pixel = pSrc->GetPixel(x, y);
      std::vector<double> spd = RGBtoSpectrum(&newColorData, *pSourceData++, *pSourceData++, *pSourceData++);
      cv::Mat opponentColor = applyTwoStageModel(&newLMSSens, spd, I, params);

      *pDestinationData++ = opponentColor.at<double>(0, 0);
      *pDestinationData++ = opponentColor.at<double>(1, 0);
      *pDestinationData++ = opponentColor.at<double>(2, 0);
    }
  }

  pDst->Convert(TMO_Yxy);
  pSrc->Convert(TMO_Yxy);

  // luminance reduction
  double epsilon = 1e-6;
  double sumLogY = 0.0;
  int pixelCount = pSrc->GetHeight() * pSrc->GetWidth();
  double Ymax = 0.0;

  // compute sum of logarithms of luminance and maximum luminance
  for (int y = 0; y < pSrc->GetHeight(); y++)
  {
    for (int x = 0; x < pSrc->GetWidth(); x++)
    {
      double Y = pSrc->GetPixel(x, y)[0];
      sumLogY += std::log(Y + epsilon);
      if (Y > Ymax)
      {
        Ymax = Y;
      }
    }
  }

  // compute average luminance
  double YLogAvg = std::exp(sumLogY / pixelCount);

  // go through the image and apply luminance reduction
  for (int y = 0; y < pDst->GetHeight(); y++)
  {
    for (int x = 0; x < pDst->GetWidth(); x++)
    {
      double Y = pSrc->GetPixel(x, y)[0];
      double Yr = luminanceReduction(Y, YLogAvg, Ymax);
      pDst->GetPixel(x, y)[0] = Yr;
    }
  }

  pDst->Convert(TMO_RGB);

  for (int y = 0; y < pDst->GetHeight(); y++)
  {
    for (int x = 0; x < pDst->GetWidth(); x++)
    {
      double *pixel = pDst->GetPixel(x, y);
      pixel[0] = pixel[0] * rcf;
      pixel[1] = pixel[1] * rcf;
      pixel[2] = pixel[2] * rcf;
    }
  }

  return 0;
}
