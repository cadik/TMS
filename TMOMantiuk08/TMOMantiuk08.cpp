/* --------------------------------------------------------------------------- *
 * TMOMantiuk08.cpp: implementation of the TMOMantiuk08 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOMantiuk08.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMantiuk08::TMOMantiuk08()
{
	SetName(L"Mantiuk08");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

#define vm 7
#define vm_contrast 2
#define round_int( x ) (int)((x) + 0.5)
#define MIN_VALUE 1e-8f
#define MAX_VALUE 1e8f
#define rgb2luminance(R,G,B) (R*0.212656 + G*0.715158 + B*0.072186)
int log_lum_cnt = round_int((8.f+8.f)/0.1f) + 1;

typedef vector< vector<double> > PixelDoubleMatrix;

struct DisplaySize{
   float view_dist;
   float pp_dist;
};

struct DisplayFunc{
   float gamma;
   float L_black;
   float L_max;
   float L_refl;
   float k;
   float E_amb;

};
struct CPDfunction{
   double min;
   double max;
   double g_max;
   double delta;
   int x_cnt;
   int g_cnt;
   int freq_cnt;
   double final_value;
   vector<double> C_val;
   vector<double> log_lum_scale;
   vector<double> g_scale;
   vector<double> f_scale;
};
struct CSFvals{
   double delta;
   vector<double> y_i;
};
typedef std::vector<CSFvals> CSFvalues;
//-------------------------------------------CDP initialization------------------------------------------------------------
void CPDinit(CPDfunction& C)
{
   C.delta = 0.1f;
   C.max = 8.f;
   C.min = -8.f;
   C.g_max = 0.7f;
   C.x_cnt = log_lum_cnt;
   C.g_cnt = round_int(C.g_max/C.delta)*2 + 1;

   int tmp;
   for(tmp = 0; tmp < 8; tmp++)
   {
      float tmp_freq = 0.5f * 30.f/(float)(1 * pow(2,tmp));
      if(tmp_freq <= 3)
      {
         break;
      }
   }
   C.freq_cnt = tmp + 1;
   for(int i=0; i < C.x_cnt*C.freq_cnt*C.g_cnt; i++)
   {
      C.C_val.push_back(0.0);
   }
   for(int i=0; i < C.x_cnt; i++)
   {
      C.log_lum_scale.push_back(C.min + C.delta*i);
   }
   for(int i=0; i < C.g_cnt; i++)
   {
      C.g_scale.push_back(-C.g_max + C.delta*i);
   }
   for(int i=0; i < C.freq_cnt; i++)
   {
      C.f_scale.push_back(0.5f * 30.f/(float)(1 * pow(2,i)));
   }
}

double calcCval(int x, int y, int z, CPDfunction& C)
{
   if((x + y*C.x_cnt + z*C.g_cnt >= 0) && (x + y*C.x_cnt + z*C.x_cnt*C.g_cnt < C.x_cnt*C.g_cnt*C.freq_cnt))
   {
      return (x + C.x_cnt + C.freq_cnt*C.g_cnt);
   }
}

//---------------------------------------------Display size functions----------------------------------------------------------------------

void displaySize(int res, float vd_screen, float vd_size, DisplaySize& disp)
{
   disp.view_dist = vd_size;
   disp.pp_dist = res * M_PI / (360 * atan(0.5/vd_screen));
}

//---------------------------------------------Display function functions-------------------------------------------------------------------
void DisplayFuncInit(float gamma, float L_black, float L_max, float k, float E_amb, DisplayFunc& df)
{
   df.gamma = gamma;
   df.L_black = L_black;
   df.L_max = L_max;
   df.k = k;
   df.E_amb = E_amb;
   df.L_refl = k/M_PI * E_amb;
}
float calcDisplayFunc(float p, DisplayFunc& df)
{
   if(p >= 0 && p <= 1)
   {
      return pow(p, df.gamma) * (df.L_max - df.L_black) + df.L_black + df.L_refl;
   }
}
float calcInverseDisplayFunc(float p, DisplayFunc& df)
{
   if(p < df.L_refl){ p = df.L_refl;}
   if(p > df.L_refl + df.L_max){ p = df.L_refl + df.L_max;}
   return pow((p - df.L_refl)/(df.L_max - df.L_black), 1/df.gamma);

}

//---------------------------------------------Human visual system functions----------------------------------------------------------------
double transducer(double G, double sens)
{
   if(vm & vm_contrast)
   {
      double W = pow(10, fabs(G)) - 1;
      double k = 0.2599, q = 3, a = 3.291, b = 3.433, e = 0.8;
      double SW = W * sens;
      int sign = 0;
      if(G < 0)
      {
         sign = -1;
      }
      else{
         sign = 1;
      }
      return sign * (a*(pow(1+pow(SW,q),1.0/3.0) - 1))/(k * pow(b + SW, e));
   }
}

//---------------------------------------------Daly's contrast sensitivity function--------------------------------------------------------
double cs_daly(double rho, double img_size, double theta, double adapt_lum, double view_dist = 0.5)
{
   double P = 250.f;
   double eps = 0.9;
   double A = 0.801*(pow(1.0 + 0.7 * 1/adapt_lum, -0.2));
   double B = 0.3 * (pow(1.0 + 100 * 1/adapt_lum, 0.15));
   double r_a = 0.856 * powf(view_dist, 0.14);
   double c = 0.0;
   double r_c = 1.0/(1.0 + 0.24 * c);
   double r_theta = 0.11 * cosf(4.0 * theta) + 0.89;
   double b_eps_rho = B * eps * rho;
   double S1 = pow(pow(3.23 * pow(rho*rho*img_size,-0.3),5)+1.0, -0.2) * A * eps * rho * exp(-b_eps_rho) * sqrt(1 + 0.06*exp(b_eps_rho));
   double new_rho = rho / (r_a * r_c * r_theta);
   double b_eps_newrho = B * eps * new_rho;
   double S2 = pow(pow(3.23 * pow(new_rho*new_rho*img_size,-0.3),5)+1.0, -0.2) * A * eps * new_rho * exp(-b_eps_newrho) * sqrt(1 + 0.06*exp(b_eps_newrho));
   if(S1 > S2)
   {
      return S2 * P;
   }
   else{
      return S1 * P;
   }
}
//-----------------------------------------------Function for computing gaussian levels-----------------------------------------------------
void GaussianLevel(int width, int height, PixelDoubleMatrix& in, PixelDoubleMatrix& out, PixelDoubleMatrix& tmp_matrix, int kernel_level){
   float a = 0.4;
   int kernel_size = 5;
   int half_kernel_size = kernel_size/2;
   float kernels[kernel_size] = {0.25 - a/2, 0.25, a, 0.25, 0.25 - a/2};

   int step = 1 * pow(2,kernel_level);

   //rows filter
   for(int i=0; i < height; i++)
   {
      for(int j=0; j < width; j++)
      {
         float tmp = 0.0;
         for(int k=0; k < kernel_size; k++)
         {
            int n = (k - half_kernel_size) * step + j;
            if(n < 0){
               n = -n;
            }
            if(n >= width){
               n = 2 * width - 2 - n;
            }
            tmp += in[i][n] * kernels[k];
         }
         tmp_matrix[i][j] = tmp;
      }
   }
   //cols filter
   for(int j=0; j < width; j++)
   {
      for(int i = 0; i < height; i++)
      {
         float tmp = 0.0;
         for(int k=0; k < kernel_size; k++)
         {
            int n = (k - half_kernel_size) * step + i;
            if(n < 0){
               n = -n;
            }
            if(n >= height){
               n = 2 * height - 2 -n;
            }
            tmp += tmp_matrix[n][j] * kernels[k];
         }
         out[i][j] = tmp;
      }
   }

};
TMOMantiuk08::~TMOMantiuk08()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOMantiuk08::Transform()
{
   // Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	pDst->Convert(TMO_Yxy); // x, y as color information
   int imageHeight = pSrc->GetHeight();
   int imageWidth = pSrc->GetWidth();
	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components

   //---------------------------------------Declarations of variables---------------------------------------------------------------------
   DisplaySize ds;
   DisplayFunc df;
   CPDfunction C;
   float res = 1024;
   float vd_screen = 2;
   float vd_size = 0.5;
   displaySize(res, vd_screen, vd_size, ds);
   DisplayFuncInit(2.2f, 0.8, 200, 0.01, 60, df);
   CPDinit(C);
   double *imgSrcData = pSrc->GetData();
   double pixelR, pixelG, pixelB;
   float threshold = 0.0043;
   PixelDoubleMatrix LogLuminancePixels(imageHeight, vector<double>(imageWidth, 0.0));
   PixelDoubleMatrix gausianOutVal(imageHeight, vector<double>(imageWidth, 0.0));
   PixelDoubleMatrix tmp_matrix(imageHeight, vector<double>(imageWidth, 0.0));
   //----------------------------------------Log luminances of pixels----------------------------------------------------------------------
   for(int j=0; j < imageHeight; j++)
   {
      for(int i=0; i<imageWidth; i++)
      {
         pixelR = *imgSrcData++;
         pixelG = *imgSrcData++;
         pixelB = *imgSrcData++;
         LogLuminancePixels[j][i] = log10(rgb2luminance(pixelR, pixelG, pixelB));
      }
   }
   bool check_out = false;
   C.final_value = 0;


   //-----------------------------------------Conditional density calculation---------------------------------------------------------------
   for(int i=0; i < C.freq_cnt; i++)
   {
      GaussianLevel(imageWidth, imageHeight, LogLuminancePixels, gausianOutVal, tmp_matrix, i);
      int g_tp = C.g_cnt/2+1;
      int g_tn = C.g_cnt/2-1;
      int g_t = C.g_cnt/2;
      for(int p=0; p < imageHeight; p++)
      {
         for(int n=0; n < imageWidth; n++)
         {
            float g = LogLuminancePixels[p][n] - gausianOutVal[p][n];
            int xi = round_int((gausianOutVal[p][n] - C.min)/C.delta);
            int gi = round_int((g + C.g_max)/C.delta);
            if(gi < 0 || gi > C.g_cnt)
            {
               continue;
            }
            if(g > threshold && g < C.delta/2)
            {
               C.C_val[calcCval(xi,g_tp,i,C)] = C.C_val[calcCval(xi,g_tp,i,C)] + 1.0;
            }
            else if(g < -threshold && g > -C.delta/2){
               C.C_val[calcCval(xi,g_tn,i,C)] = C.C_val[calcCval(xi,g_tn,i,C)] + 1.0;
            }
            else{
               C.C_val[calcCval(xi,gi,i,C)] = C.C_val[calcCval(xi,gi,i,C)] + 1.0;
            }
         }
      }
      for(int m=0; m < C.x_cnt; m++)
      {
         if(C.C_val[calcCval(m,g_t,i,C)] == 0)
         {
            continue;
         }
         bool grad = false;
         for(int n=0; n < C.g_cnt; n++)
         {
            if(n != g_t && C.C_val[calcCval(m,n,i,C)] != 0)
            {
               grad = true;
               break;
            }
         }
         if(!grad)
         {
            C.C_val[calcCval(m,g_tp,i,C)] = C.C_val[calcCval(m,g_tp,i,C)] + 1.0;
            C.C_val[calcCval(m,g_tn,i,C)] = C.C_val[calcCval(m,g_tn,i,C)] + 1.0;
         }
         for(int k=0; k < C.g_cnt; k++)
         {
            C.final_value += C.C_val[calcCval(m,k,i,C)];
         }
      }
   }
   fprintf(stderr,"final val %g\n",C.final_value);
   fprintf(stderr,"size %d\n",C.C_val.size());


   
	//------------------------------------------Solving quadratic problem-----------------------------------------------
   double display_drange = log10(calcDisplayFunc(1.0f, df)/calcDisplayFunc(0.0f, df));
   CSFvalues csf;
   for(int f=0; f< C.freq_cnt; f++)
   {
      CSFvals tmp;
      for(int i=0; i < C.x_cnt; i++)
      {
         tmp.y_i.push_back(cs_daly(C.f_scale[f],1,0,pow(10.0,C.log_lum_scale[i])));
         tmp.delta = C.log_lum_scale[1] - C.log_lum_scale[0];
      }
      csf.push_back(tmp);
   }














	double pY, px, py;

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
	return 0;
}
