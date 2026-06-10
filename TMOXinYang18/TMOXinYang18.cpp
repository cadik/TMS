/* --------------------------------------------------------------------------- *
 * TMOXinYang18.cpp: implementation of the TMOXinYang18 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOXinYang18.h"
#include <tensorflow/c/c_api.h>
//#include <Python.h>
#include <iostream>
//moze byt na pozadi spusteny server?
//#include <triton/core/tritonserver.h>
//#include <triton/core/tritonserver_common.h>
//#include <triton/backend/backend_common.h>


using namespace std;
/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOXinYang18::TMOXinYang18()
{
	SetName(L"XinYang18");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOXinYang18::~TMOXinYang18()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
//* int TMOXinYang18::Transform()
//{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
//*	pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
//	pDst->Convert(TMO_Yxy); // x, y as color information

//	double *pSourceData = pSrc->GetData();		// You can work at low level data
//	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
/*												// three colour components
	double pY, px, py;
   cerr<< "fafaf"<< endl;
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			// Here you can use your transform
			// expressions and techniques...
			pY *= dParameter; // Parameters can be used like
							  // simple variables

			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);

   //./tmocmd XinYang18 ~/Desktop/bakalarka/TMS/TMOXinYang18/DRHT-python/input/gau_IMG_3238.JPG

   //tl::Tensor x;
/*

   FILE* file;
    int argc;
    wchar_t * argv[4];

   /* argc = 4;
    argv[0] = L"/home/jalinux/Desktop/bakalarka/TMS/TMOXinYang18/DRHT-python/main.py";
    argv[1] = L"/home/jalinux/Desktop/bakalarka/TMS/TMOXinYang18/DRHT-python/input/gau_IMG_3238.JPG";
    argv[2] = L"2448";
    argv[3] = L"3264";
    Py_SetPythonHome(L"/usr/lib/python3.8");
    Py_SetProgramName(argv[0]);
    Py_Initialize();
    PySys_SetArgv(argc, argv);
    file = fopen("/home/jalinux/Desktop/bakalarka/TMS/TMOXinYang18/DRHT-python/main.py","r");
    cerr<<file<<" <- this"<<endl;
    PyRun_SimpleFile(file, "/home/jalinux/Desktop/bakalarka/TMS/TMOXinYang18/DRHT-python/main.py");
    cerr<<"HRR";
    Py_Finalize();

    fclose(file);

   //PyObject* ModuleString = PyUnicode_FromString("~/Desktop/bakalarka/TMS/TMOXinYang18/DRHT-python/img_io.py");
   //PyObject *pModule = PyImport_Import(ModuleString);
   //PyObject *pModule = PyImport_ImportModule("~/Desktop/bakalarka/TMS/TMOXinYang18/DRHT-python/img_io.py");//////////////////////
   

   

	return 0;
}

tl::Tensor model(const tl::Tensor& x, const tl::Tensor& weights) {
    // Convert input tensor values to [0, 255] range
    tl::Tensor x_in = tl::Multiply(x, tl::Scalar(255.0));

    // Create an InputLayer
    tl::InputLayer net_in(x_in, "input_layer");

    // Run the encoder
    std::vector<tl::Layer*> conv_layers, skip_layers;
    std::tie(conv_layers, skip_layers) = encoder(net_in, weights);

    // Create the Conv2dLayer
    tl::Conv2dLayer network(conv_layers,
                             other arguments ,
                            "encoder/h6/conv");

    // Create the BatchNormLayer
    tl::BatchNormLayer network_batch_norm(network, false, "encoder/h6/batch_norm");

    // Apply ReLU activation
    tl::Scalar zero(0.0);
    network_batch_norm.outputs = tl::ReLU(network_batch_norm.outputs, zero, "encoder/h6/relu");

    // Run the decoder
    network = decoder(network_batch_norm, skip_layers);

    return network;
}

class Layer{
   public:
      string name;
      Layer(string name = nullptr){
         if()
      }
};

class InputLayer : public Layer{

};

void model(float* x,int cnt){
      for(int i = 0; i<cnt;i++) x[i]*=255;
      
      //float net_in = 

      //float network = decoder()
}*/
