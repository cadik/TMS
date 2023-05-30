/* --------------------------------------------------------------------------- *
 * TMOJah06.cpp: implementation of the TMOJah06 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOJah06.h"

 /* --------------------------------------------------------------------------- *
  * Constructor serves for describing a technique and input parameters          *
  * --------------------------------------------------------------------------- */
TMOJah06::TMOJah06()
{
   SetName(L"Jah06");
   SetDescription(L"apply style transfer using neural network");

   /*
   dPath.SetName(L"p");
   dPath.SetDescription(L"Path to torchscript model");
   dPath.SetDefault("model.pt");
   this->Register(dPath);
   */
   //cv::useOptimized();
}

TMOJah06::~TMOJah06()
{
}

cv::Mat TMOImage2Mat(TMOImage* pSrc)
{
   double* pSourceData;
   int rowsCnt, colsCnt;

   pSourceData = pSrc->GetData();
   rowsCnt = pSrc->GetHeight();
   colsCnt = pSrc->GetWidth();

   cv::Mat srcConvMat(rowsCnt, colsCnt, CV_64FC3);

   if (srcConvMat.isContinuous())
   {
      colsCnt *= rowsCnt;
      rowsCnt = 1;
   }
   int chCnt = 3;
   for (int y = 0; y<rowsCnt; y++)
   {
      for (int x = 0; x<colsCnt; x++)
      {
         for (int c = 0; c<chCnt; c++)
         {
            srcConvMat.at<double>(y, x*chCnt+c) = pSourceData[chCnt-c-1];
         }

         /* Add count of channels (RGB) to pointer */
         pSourceData += chCnt;
      }
   }

   return srcConvMat;
}

// Function to load the PyTorch model
torch::jit::script::Module loadModel(const std::string& modelPath) {
   torch::jit::script::Module module;
   try {
      // Load the model
      module = torch::jit::load(modelPath);
   }
   catch (const c10::Error& e) {
      std::cerr<<"Error loading the model: "<<e.what()<<std::endl;
      exit(1);
   }
   return module;
}

// Function to preprocess the input image
torch::Tensor preprocessImage(const cv::Mat& img) {
   cv::Mat floatImage;
   img.convertTo(floatImage, CV_32FC3, 1.0/255.0);

   //convert to tensor
   torch::Tensor tensor = torch::from_blob(floatImage.data, { 1, floatImage.rows, floatImage.cols, 3 }, torch::kFloat32);
   tensor = tensor.permute({ 0, 3, 1, 2 }); // Transpose dimensions from [1, H, W, C] to [1, C, H, W]

   return tensor;
}

cv::Mat tensorToImage(const torch::Tensor& tensor) {
   // Remove the batch_size dimension
   torch::Tensor squeezed = tensor.squeeze();

   // Convert tensor to CPU and extract data as a pointer
   torch::Tensor cpuTensor = squeezed.to(torch::kCPU);
   float* dataPtr = cpuTensor.data_ptr<float>();

   // Get the dimensions of the tensor
   int64_t channels = cpuTensor.size(0);
   int64_t height = cpuTensor.size(1);
   int64_t width = cpuTensor.size(2);

   // Create an OpenCV Mat with the same dimensions
   cv::Mat image(height, width, CV_32FC(channels));

   // Copy the data from the tensor to the Mat
   std::memcpy(image.data, dataPtr, channels*height*width*sizeof(float));

   // Transpose from [C, H, W] -> [H, W, C]
   cv::Mat transposedImage = image.t();

   // Convert the image to the appropriate data type and scale
   cv::Mat scaledImage;
   transposedImage.convertTo(scaledImage, CV_8UC(channels));

   // Return the final image
   return scaledImage;
}
/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOJah06::Transform()
{
   // model loading
   torch::jit::script::Module model = loadModel("model.pt");

   // image load + preprocess
   cv::Mat inputImage = TMOImage2Mat(pSrc);
   torch::Tensor inputTensor = preprocessImage(inputImage);
   inputTensor = inputTensor.to(torch::kCPU);

   torch::Tensor generatedTensor = model.forward({ inputTensor }).toTensor();
   cv::Mat generatedImage = tensorToImage(generatedTensor);


   //image save
   double* pDestinationData = pDst->GetData();

   int chCnt = 3;
   int y = 0;
   for (y = 0; y<pSrc->GetHeight(); y++)
   {
      pSrc->ProgressBar(y, pSrc->GetHeight());
      for (int x = 0; x<pSrc->GetWidth(); x++)
      {
         // Convert the pixel values to the appropriate data type
         for (int c = 0; c<chCnt; c++) {
            *pDestinationData++ = static_cast<double>(generatedImage.at<cv::Vec3f>(y, x)[chCnt-c-1]);
         }
      }
   }
   pSrc->ProgressBar(y, pSrc->GetHeight());
   return 0;
}
