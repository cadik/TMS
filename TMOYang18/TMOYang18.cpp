/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2023-24                                           *
*                                                                              *
*                       Implementation of the TMOYang18 class            *
*                                                                              *
*                       Author: Lukas Macejka (xmacej03)                       *
*                       Mail: xmacej03@vutbr.cz                                *
*                                                                              *
*******************************************************************************/
#include "TMOYang18.h"
#include <opencv2/opencv.hpp>
#include <tensorflow/c/c_api.h>
#include <iostream>
#include <set>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>

using namespace std;
using namespace cv;

/**
 * @brief Constructor serves for describing a technique and input parameters
 */
TMOYang18::TMOYang18()
{
    SetName(L"Yang18");
    SetDescription(L"Image Correction via Deep Reciprocating HDR Transformation");

    dParameter.SetName(L"ParameterName");
    dParameter.SetDescription(L"ParameterDescription");
    dParameter.SetDefault(1);
    dParameter = 1.;
    dParameter.SetRange(-1000.0, 1000.0);
    this->Register(dParameter);
}

/**
 * @brief Function to read the contents of a file into a TensorFlow buffer (loading binary data)
 * 
 * @param file_path Path to the file to be read
 * @param status TensorFlow status object for error handling
 * @return TF_Buffer* Pointer to the TensorFlow buffer containing the file data
 */
TF_Buffer* read_file_into_tf_buffer(const char* file_path, TF_Status* status) {
    std::ifstream file(file_path, std::ios::binary | std::ios::ate);
    if (!file) {
        TF_SetStatus(status, TF_INVALID_ARGUMENT, "Nejde otvorit subor");
        return nullptr;
    }

    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size)) {
        TF_SetStatus(status, TF_INVALID_ARGUMENT, "Nejde nacitat subor");
        return nullptr;
    }

    return TF_NewBufferFromString(buffer.data(), buffer.size());
}

/**
 * @brief Function to execute a TensorFlow graph (model)
 * 
 * @param input_tensor Input tensor to be processed
 * @param graph TensorFlow graph containing the model
 * @param status TensorFlow status object for error handling
 * @param session_opts TensorFlow session options
 * @param session TensorFlow session
 * @return float* Pointer to the output data from the model
 */
float* ExecuteCNN(TF_Tensor* input_tensor, TF_Graph* graph, TF_Status* status, TF_SessionOptions* session_opts, TF_Session* session) {
    TF_Operation* input_ops = TF_GraphOperationByName(graph, "Placeholder");
    TF_Operation* output_ops = TF_GraphOperationByName(graph, "add_7");
    if (input_ops == nullptr || output_ops == nullptr) {
        fprintf(stderr, "Error: Unable to find input or output operation\n");
    }

    TF_Output input_op = {input_ops, 0};
    TF_Output output_op = {output_ops, 0};
    const TF_Output inputs[] = {input_op};
    TF_Tensor* outputs[1] = {nullptr}; 
    
    TF_SessionRun(session, nullptr, inputs, &input_tensor, 1, &output_op, outputs, 1, nullptr, 0, nullptr, status);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "Error running inference: %s\n", TF_Message(status));
    }

    float* output_data = static_cast<float*>(TF_TensorData(outputs[0]));
    
    return output_data;
}

/**
 * @brief Helper function to resize an image 
 * 
 * @param inputImage Input image to be resized
 * @param newWidth New width for the resized image
 * @param newHeight New height for the resized image
 * @return Mat Resized image
 */
Mat resizeImage(Mat inputImage, int newWidth, int newHeight) {
    Mat resizedImage;
    resize(inputImage, resizedImage, Size(newWidth, newHeight), 0, 0, INTER_LINEAR);
    return resizedImage;
}

/**
 * @brief Converts TMOImage(input) to cv::Mat
 * 
 * @param width Width of the image
 * @param height Height of the image
 * @param input Input image in TMOImage format
 * @return Mat Image in cv::Mat format
 */
Mat TMOImage2Mat(int width, int height, Mat *input) {
    int rowsCnt, colsCnt;

    rowsCnt = height;
    colsCnt = width;

    Mat srcConvMat(rowsCnt, colsCnt, CV_32FC3);
    
    Mat tempMat(rowsCnt*colsCnt,1,CV_32FC3);
    int i = 0;
    
    for (int x = 0; x < colsCnt; x++) {
        for (int y = 0; y < rowsCnt; y++) {
            double R = input->at<float>(0,y*colsCnt + x);
            double G = input->at<float>(1,y*colsCnt + x);
            double B = input->at<float>(2,y*colsCnt + x);
            tempMat.at<Vec3f>(i,0) = Vec3f(B,G,R); // BGR!!!
            i++;
        }
    }

    i=0;
    for (int y = 0; y < rowsCnt ; y++) {
        for (int x = 0; x <  colsCnt ; x++) {
            srcConvMat.at<Vec3f>(y,x) = tempMat.at<Vec3f>(i,0);
            i++;
        }
    }
    return srcConvMat;
}

/**
 * @brief Converts back cv::Mat to TMOImage for output
 * 
 * @param width Width of the image
 * @param height Height of the image
 * @param input Input image in cv::Mat format
 * @return Mat Image in TMOImage format
 */
Mat Mat2TMOImage(int width, int height,  Mat *input) {
    int rowsCnt = height;
    int colsCnt = width;

    Mat out(3, rowsCnt * colsCnt, CV_32FC3);

    Mat outh(width*height, 1, CV_32FC3);

    unsigned i=0;
    for (int y = 0; y < rowsCnt ; y++) {
        for (int x = 0; x <  colsCnt ; x++) {
            outh.at<Vec3f>(i,0) = input->at<Vec3f>(y,x);
            i++;
        }
    }
    i=0;
    for (int x = 0; x < colsCnt; x++) {
        for (int y = 0; y < rowsCnt; y++) {
            double r = outh.at<Vec3f>(i,0)[0];
            double g = outh.at<Vec3f>(i,0)[1];
            double b = outh.at<Vec3f>(i,0)[2];
            out.at<float>(0,y*colsCnt + x) = b;
            out.at<float>(1,y*colsCnt + x) = g;
            out.at<float>(2,y*colsCnt + x) = r;       
            i++;
        }
    }
    return out;
}

/**
 * @brief Destructor for the TMOYang18 class
 */
TMOYang18::~TMOYang18()
{
}

/**
 * @brief This overloaded function is an implementation of the tone mapping operator
 * 
 * @return int Status code
 */
int TMOYang18::Transform()
{
    //std::cerr<< "Hello from TensorFlow C library version "<< TF_Version() << std::endl;

    TF_Graph* graph = TF_NewGraph();
    TF_Status* status = TF_NewStatus();

    TF_Buffer* graph_def = read_file_into_tf_buffer("../TMOYang18/DRHT-python/checkpoint/hdr2ldr/model.pb", status);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "Chyba nacitania grafu: %s\n", TF_Message(status));
    }

    TF_ImportGraphDefOptions* import_opts = TF_NewImportGraphDefOptions();
    TF_GraphImportGraphDef(graph, graph_def, import_opts, status);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "Chyba pri importovani do grafu: %s\n", TF_Message(status));
    }

    TF_DeleteImportGraphDefOptions(import_opts);
    TF_DeleteBuffer(graph_def);

    TF_SessionOptions* session_opts = TF_NewSessionOptions();
    TF_Session* session = TF_NewSession(graph, session_opts, status);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "Chyba pri vytvarani session: %s\n", TF_Message(status));
    }

    double *pSourceData = pSrc->GetData();        
    double *pDestinationData = pDst->GetData(); 
    float R, G, B;
    
    int pixels_c = pSrc->GetHeight() * pSrc->GetWidth();
    Mat RGB_M(3,pixels_c,CV_32FC1);
    double *p_source_tmp = pSourceData;

    for(int i=0;i<pSrc->GetWidth();i++){
        for(int j=0;j<pSrc->GetHeight();j++){
            for(int k=0;k<3;k++){
                RGB_M.at<float>(k,(j*pSrc->GetWidth())+i) = (*p_source_tmp);
                p_source_tmp++;
            }
        }
    }

    std::vector<float> img_in_array;

    Mat tester = TMOImage2Mat(pSrc->GetWidth(), pSrc->GetHeight(),&RGB_M);

    Mat resizedImage = resizeImage(tester, 1024, 768);
    
    Mat transposedImage = Mat2TMOImage(1024, 768, &resizedImage);

    int rows = transposedImage.rows;
    int cols = transposedImage.cols;

    for(int i=0;i<1024;i++){
        for(int j=0;j<768;j++){
            B = transposedImage.at<float>(2,(j*1024)+i);
            G = transposedImage.at<float>(1,(j*1024)+i);
            R = transposedImage.at<float>(0,(j*1024)+i);
            img_in_array.push_back(R);
            img_in_array.push_back(G);
            img_in_array.push_back(B);
        }
    }

    const std::vector<int64_t> dims = {1, 768, 1024, 3};
    TF_Tensor* input_tensor = TF_NewTensor(TF_FLOAT, dims.data(), dims.size(), img_in_array.data(),
                                            img_in_array.size()*sizeof(float), [](void* data, size_t, void* arg) {}, nullptr);

    float* returning_data = ExecuteCNN(input_tensor, graph, status, session_opts, session);
    int j,f = 0;
    Mat ProPixelOut(3,1024*768,CV_32FC1);

    for (int i = 0; i < 1024; i++) {
        for (j = 0; j < 768; j++) {
            ProPixelOut.at<float>(0,(j*1024)+i) = returning_data[f];
            f++;
            ProPixelOut.at<float>(1,(j*1024)+i) = returning_data[f];
            f++;
            ProPixelOut.at<float>(2,(j*1024)+i) = returning_data[f];
            f++;
        }
    }

    Mat CNNOut1 = TMOImage2Mat(1024, 768,&ProPixelOut);

    Mat CNNOut2 = resizeImage(CNNOut1, pSrc->GetWidth(), pSrc->GetHeight());

    Mat CNNOut3 = Mat2TMOImage(pSrc->GetWidth(), pSrc->GetHeight(), &CNNOut2);

    for(int i=0;i<pSrc->GetWidth();i++){
        for(int j=0;j<pSrc->GetHeight();j++){
            for(int k=0;k<3;k++){
                *pDestinationData++ =  CNNOut3.at<float>(k,(j*pSrc->GetWidth())+i);
            }
        }
    }

    pSrc->ProgressBar(j, 768);
    pDst->Convert(TMO_RGB);

    TF_DeleteSession(session, status);
    TF_DeleteSessionOptions(session_opts);
    TF_DeleteGraph(graph);
    TF_DeleteStatus(status);
    return 0;
}
