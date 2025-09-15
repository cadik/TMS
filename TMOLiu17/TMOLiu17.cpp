/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*          Log-Euclidean Metrics for Contrast Preserving Decolorization        *
*                                                                              *
*             Author: Peter Zdravecký [xzdrav00 AT stud.fit.vutbr.cz]          *
*                                    Brno 2024                                 *
*                                                                              *
*******************************************************************************/

#include "TMOLiu17.h"

TMOLiu17::TMOLiu17() {
  SetName(L"Liu17");
  SetDescription(L"Log-Euclidean Metrics for Contrast Preserving Decolorization");

  gamma.SetName(L"g");
  gamma.SetDescription(L"Gamma value");
  gamma.SetDefault(0.2);
  gamma = 0.2;
  this->Register(gamma);
}

TMOLiu17::~TMOLiu17() {}

/**
 * @brief Generates a random permutation of indices.
 *
 * @param n The number of indices.
 * @param indices The output matrix to store the shuffled indices.
 */
void randperm(int n, Mat &indices) {
  indices = Mat::zeros(n, 1, CV_32S);
  for (int i = 0; i < n; ++i) {
    indices.at<int>(i, 0) = i;
  }
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(indices.begin<int>(), indices.end<int>(), g);
}

int TMOLiu17::Transform() {
  pSrc->Convert(TMO_RGB);
  pDst->Convert(TMO_RGB);
  double *pSourceData = pSrc->GetData();
  double *pDestinationData = pDst->GetData();

  double gammaValue = gamma;
  double Lpp = 0.013;

  int height = pSrc->GetHeight();
  int width = pSrc->GetWidth();

  Mat im = Mat(height, width, CV_64FC3, pSourceData);
  Mat W = wei();

  int n = im.rows;
  int m = im.cols;
  Size newShape(round(64 / sqrt(n * m) * n), round(64 / sqrt(n * m) * m));

  Mat ims;
  resize(im, ims, newShape, 0, 0, INTER_NEAREST_EXACT);

  Mat R, G, B;
  Mat channels[3];
  split(ims, channels);
  R = channels[0];
  G = channels[1];
  B = channels[2];

  Mat imV;
  hconcat(R.reshape(1, R.total()), G.reshape(1, G.total()), imV);
  hconcat(imV, B.reshape(1, B.total()), imV);

  Mat t1;
  randperm(imV.rows, t1);
  Mat imV_shuffled = Mat::zeros(imV.size(), imV.type());
  for (int i = 0; i < imV.rows; ++i) {
    imV.row(t1.at<int>(i, 0)).copyTo(imV_shuffled.row(i));
  }
  Mat Pg = imV - imV_shuffled;

  resize(ims, ims, Size(ims.cols / 2, ims.rows / 2), 0, 0, INTER_NEAREST_EXACT);

  Mat Rx, Gx, Bx;
  Mat Ry, Gy, By;
  Mat channels2[3];
  split(ims, channels2);
  Rx = channels2[0].colRange(0, channels2[0].cols - 1) -
       channels2[0].colRange(1, channels2[0].cols);
  Gx = channels2[1].colRange(0, channels2[1].cols - 1) -
       channels2[1].colRange(1, channels2[1].cols);
  Bx = channels2[2].colRange(0, channels2[2].cols - 1) -
       channels2[2].colRange(1, channels2[2].cols);
  Ry = channels2[0].rowRange(0, channels2[0].rows - 1) -
       channels2[0].rowRange(1, channels2[0].rows);
  Gy = channels2[1].rowRange(0, channels2[1].rows - 1) -
       channels2[1].rowRange(1, channels2[1].rows);
  By = channels2[2].rowRange(0, channels2[2].rows - 1) -
       channels2[2].rowRange(1, channels2[2].rows);

  Mat RxGxBx;
  hconcat(Rx.reshape(1, Rx.total()), Gx.reshape(1, Gx.total()), RxGxBx);
  hconcat(RxGxBx, Bx.reshape(1, Bx.total()), RxGxBx);
  Mat RyGyBy;
  hconcat(Ry.reshape(1, Ry.total()), Gy.reshape(1, Gy.total()), RyGyBy);
  hconcat(RyGyBy, By.reshape(1, By.total()), RyGyBy);
  Mat Pl;
  vconcat(RxGxBx, RyGyBy, Pl);

  Mat PT;
  vconcat(Pg, Pl, PT);

  Mat det;
  reduce(PT.mul(PT), det, 1, REDUCE_SUM);
  sqrt(det, det);
  det /= 1.41;

  Mat mask = det < 0.05;
  Mat filtered_PT, filtered_det;
  PT.copyTo(filtered_PT, ~repeat(mask, 1, 3));
  det.copyTo(filtered_det, ~mask);

  Mat detM;
  repeat(filtered_det.reshape(1, filtered_det.total()), 1, W.rows, detM);

  Mat L = filtered_PT * W.t();

  Mat log_val = abs(L) / (abs(detM) + Lpp);
  Mat coef;
  log(log_val, coef);
  coef = abs(coef);

  Mat U;
  exp(-gamma * coef, U);

  Mat Es;
  reduce(U, Es, 0, REDUCE_AVG);

  Point maxLoc;
  minMaxLoc(Es, NULL, NULL, NULL, &maxLoc);
  int bw = maxLoc.x;

  Mat im_channels[3];
  split(im, im_channels);
  Mat g = W.at<double>(bw, 0) * im_channels[0] +
          W.at<double>(bw, 1) * im_channels[1] +
          W.at<double>(bw, 2) * im_channels[2];
  cv::normalize(g, g, 0, 1, NORM_MINMAX);

  int N = width * height;
  for (int i = 0; i < N; i++) {
    pDestinationData[i * 3] = g.at<double>(i / width, i % width);
    pDestinationData[i * 3 + 1] = g.at<double>(i / width, i % width);
    pDestinationData[i * 3 + 2] = g.at<double>(i / width, i % width);
  }
  return 0;
}

/**
 * @brief Generates the weight matrix for discrete optimization.
 *
 * @return The weight matrix.
 */
Mat TMOLiu17::wei() {
  Mat W = (Mat_<double>(66, 3) << 0.0, 0.0, 1.0, 0.0, 0.1, 0.9, 0.0, 0.2, 0.8,
           0.0, 0.3, 0.7, 0.0, 0.4, 0.6, 0.0, 0.5, 0.5, 0.0, 0.6, 0.4, 0.0, 0.7,
           0.3, 0.0, 0.8, 0.2, 0.0, 0.9, 0.1, 0.0, 1.0, 0.0, 0.1, 0.0, 0.9, 0.1,
           0.1, 0.8, 0.1, 0.2, 0.7, 0.1, 0.3, 0.6, 0.1, 0.4, 0.5, 0.1, 0.5, 0.4,
           0.1, 0.6, 0.3, 0.1, 0.7, 0.2, 0.1, 0.8, 0.1, 0.1, 0.9, 0.0, 0.2, 0.0,
           0.8, 0.2, 0.1, 0.7, 0.2, 0.2, 0.6, 0.2, 0.3, 0.5, 0.2, 0.4, 0.4, 0.2,
           0.5, 0.3, 0.2, 0.6, 0.2, 0.2, 0.7, 0.1, 0.2, 0.8, 0.0, 0.3, 0.0, 0.7,
           0.3, 0.1, 0.6, 0.3, 0.2, 0.5, 0.3, 0.3, 0.4, 0.3, 0.4, 0.3, 0.3, 0.5,
           0.2, 0.3, 0.6, 0.1, 0.3, 0.7, 0.0, 0.4, 0.0, 0.6, 0.4, 0.1, 0.5, 0.4,
           0.2, 0.4, 0.4, 0.3, 0.3, 0.4, 0.4, 0.2, 0.4, 0.5, 0.1, 0.4, 0.6, 0.0,
           0.5, 0.0, 0.5, 0.5, 0.1, 0.4, 0.5, 0.2, 0.3, 0.5, 0.3, 0.2, 0.5, 0.4,
           0.1, 0.5, 0.5, 0.0, 0.6, 0.0, 0.4, 0.6, 0.1, 0.3, 0.6, 0.2, 0.2, 0.6,
           0.3, 0.1, 0.6, 0.4, 0.0, 0.7, 0.0, 0.3, 0.7, 0.1, 0.2, 0.7, 0.2, 0.1,
           0.7, 0.3, 0.0, 0.8, 0.0, 0.2, 0.8, 0.1, 0.1, 0.8, 0.2, 0.0, 0.9, 0.0,
           0.1, 0.9, 0.1, 0.0, 1.0, 0.0, 0.0);
  return W;
}
