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

#define round_int( x ) (int)((x) + 0.5)
#define MIN_VALUE 1e-8f
#define MAX_VALUE 1e8f
#define rgb2luminance(R,G,B) (R*0.212656 + G*0.715158 + B*0.072186)
#define log_lum_cnt (round_int((8.f+8.f)/0.1f) + 1);

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
   vector<float> y_i;
   vector<float> x_i;
   size_t v_size;
};
struct ToneCurve{
   vector<float> y_i;
   vector<float> x_i;
};
typedef std::vector<CSFvals> CSFvalues;
//-------------------------------------------CDP initialization------------------------------------------------------------
float minvalue(float x)
{
   if(x > MIN_VALUE)
   {
      return x;
   }
   else
   {
      return MIN_VALUE;
   }
}
void CPDinit(CPDfunction& C)
{
   C.delta = 0.1f;
   C.max = 8.f;
   C.min = -8.f;
   C.g_max = 0.7f;
   C.x_cnt = log_lum_cnt;
   C.g_cnt = round_int(C.g_max/C.delta)*2 + 1;
   for(int i=0; i < C.x_cnt; i++)
   {
      C.log_lum_scale.push_back(0);
   }
   int tmp;
   for(tmp = 0; tmp < 8; tmp++)
   {
      float tmp_freq = 0.5f * 30.f/(float)(1*pow(2,tmp));
      if(tmp_freq <= 3)
      {
         break;
      }
   }
   C.freq_cnt = tmp + 1;
   for(int i=0; i < C.x_cnt*C.freq_cnt*C.g_cnt; i++)
   {
      C.C_val.push_back(0);
   }
   if(C.log_lum_scale[0] == 0)
   {
      for(int i=0; i < C.x_cnt; i++)
      {
         C.log_lum_scale[i] = C.min + C.delta*i;
      }
   }
   
   for(int i=0; i < C.g_cnt; i++)
   {
      C.g_scale.push_back(-1*C.g_max + C.delta*i);
   }
   for(int i=0; i < C.freq_cnt; i++)
   {
      C.f_scale.push_back(0.5f * 30.f/(float)(1<<i));
   }
}

double calcCval(int x, int y, int z, CPDfunction& C)
{
   assert((x + y*C.x_cnt + z*C.x_cnt*C.g_cnt >= 0) && (x + y*C.x_cnt + z*C.x_cnt*C.g_cnt < C.x_cnt*C.g_cnt*C.freq_cnt));
   return x + y*C.x_cnt + z*C.x_cnt*C.g_cnt;
}
//---------------------------------------------Interpolation func--------------------------------------------------------------------------
double interpolation(double val, CSFvals& csf)
{
   csf.delta = csf.x_i[1]-csf.x_i[0];
   double f = (val - csf.x_i[0])/csf.delta;
   size_t l = (size_t)(f);
   size_t h = (size_t)ceil(f);
   //fprintf(stderr,"f: %g l: %g h: %g\n",f,l,h);
   if(f < 0)
   {
      return csf.y_i[0];
   }
   if(h >= csf.v_size)
   {
      return csf.y_i[csf.v_size-1];
   }
   if(l == h)
   {
      return csf.y_i[l];
   }
   return csf.y_i[l] + (csf.y_i[h] - csf.y_i[l])*(f - (double)l);
}

//---------------------------------------------Display size functions----------------------------------------------------------------------

void displaySize(float res,float vd_size, DisplaySize& disp)
{
   disp.view_dist = vd_size;
   disp.pp_dist = res;
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
   assert(p >= 0 && p <= 1);
   return pow(p, df.gamma) * (df.L_max - df.L_black) + df.L_black + df.L_refl;
}
float calcInverseDisplayFunc(float p, DisplayFunc& df)
{
   if(p < df.L_refl + df.L_black){ p = df.L_refl + df.L_black;}
   if(p > df.L_refl + df.L_max){ p = df.L_refl + df.L_black + df.L_max;}
   return powf((p - df.L_refl - df.L_black)/(df.L_max - df.L_black), 1/df.gamma);

}

//---------------------------------------------Human visual system functions----------------------------------------------------------------
double transducer(double G, double sens)
{
   
   double W = pow(10, fabs(G)) - 1;
   double k = 0.2599, q = 3.0, a = 3.291, b = 3.433, e = 0.8;
   double SW = W * sens;
   int sign = 0;
   if(G < 0)
   {
      sign = -1;
   }
   else{
      sign = 1;
   }
   return sign * a*(pow(1.+pow(SW,q),1./3.) - 1.)/(k * pow(b + SW, e));
   //return G * sens;
}

//---------------------------------------------Daly's contrast sensitivity function--------------------------------------------------------
double cs_daly(double rho, double img_size, double theta, double adapt_lum, double view_dist = 0.5)
{
   //adapt_lum = 1000.;
   //rho = 4.;
   if(rho == 0)
   {
      return 0;
   }
   double P = 250.f;
   double eps = 0.9;
   double A = 0.801*pow(1 + 0.7 / adapt_lum, -0.20);
   double B = 0.3 * (pow(1 + 100 / adapt_lum, 0.15));
   double r_a = 0.856 * powf(view_dist, 0.14);
   double c = 0.0;
   double r_c = 1.0/(1.0 + 0.24 * c);
   double r_theta = 0.11 * cosf(4.0 * theta) + 0.89;
   double b_eps_rho = B * eps * rho;
   double S1 = pow(pow(3.23 * pow(rho*rho*img_size,-0.3),5.0)+1.0, -0.2) * A * eps * rho * exp(-b_eps_rho) * sqrt(1 + 0.06*exp(b_eps_rho));
   double new_rho = rho / (r_a * r_c * r_theta);
   double b_eps_newrho = B * eps * new_rho;
   double S2 = powf(pow(3.23 * pow(new_rho*new_rho*img_size,-0.3),5.0)+1.0, -0.2) * A * eps * new_rho * exp(-b_eps_newrho) * sqrt(1 + 0.06*exp(b_eps_newrho));
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

   int step = 1 * pow(2, kernel_level);

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

}
//-------------------------------function to calculate final y_i--------------------------------------------
void calculateToneCurve(ToneCurve& tc, vector<int>& unused, gsl_vector *x, int cnt, int xcnt,double Lmin, double Lmax)
{
   double alpha = 1;
   double sum = 0.0;
   for(int i=0; i < cnt; i++)
   {
      sum += gsl_vector_get(x,i);
   }
   double tmp = log10(Lmin) + alpha*(log10(Lmax) - log10(Lmin) - sum);
   tc.y_i[0] = tmp;
   double v = 0;
   for(int i =0; i < xcnt -1; i++)
   {
      if(unused[i] != -1)
      {
         int calc;
         for(calc = i+1; calc < (xcnt-1) && unused[calc] == -1; calc++);
         if(calc == (xcnt-1))
         {
            v = 0;
            tc.y_i[i] = tmp;
            tmp += gsl_vector_get(x, unused[i]);
            continue;
         }
         else
         {
            v = gsl_vector_get(x, unused[i]) / (double)(calc - i);
         }
      }
      tc.y_i[i] = tmp;
      tmp += v;
   }
   tc.y_i[xcnt-1] = tmp;
}
void multiple(gsl_matrix *a, gsl_matrix *b, gsl_vector *x)
{
   assert(a->size1 == x->size);
   for(int i=0; i < a->size2; i++)
   {
      for(int j=0; j < a->size1; j++)
      {
         gsl_matrix_set(b,j,i,gsl_matrix_get(a,j,i)*gsl_vector_get(x,j));
      }
   }
}

class auto_cqpminimizer
{
   gsl_cqpminimizer *v;
   public:
      auto_cqpminimizer(gsl_cqpminimizer *v): v(v){}
      ~auto_cqpminimizer() { gsl_cqpminimizer_free(v);}
      operator gsl_cqpminimizer* () {return v;}
};
const static gsl_matrix null_matrix = {0};
const static gsl_vector null_vector = {0};
void solver(gsl_matrix *Q, gsl_vector *q, gsl_matrix *C, gsl_vector *d, gsl_vector *x)
{
   gsl_cqp_data cqpd;
   cqpd.Q = Q;
   cqpd.q = q;
   cqpd.A = &null_matrix;
   cqpd.b = &null_vector;
   cqpd.C = C;
   cqpd.d = d;
   size_t n = cqpd.Q->size1;
   size_t me = cqpd.b->size;
   size_t mi = cqpd.d->size;
   const size_t max_iter = 100;
   size_t iter = 1;
   int status;
   const gsl_cqpminimizer_type * T;
   T = gsl_cqpminimizer_mg_pdip;
   auto_cqpminimizer s(gsl_cqpminimizer_alloc(T, n, me, mi));
   status = gsl_cqpminimizer_set(s, &cqpd);
   bool verbose = false;
   do
   {
      status = gsl_cqpminimizer_iterate(s);
      status = gsl_cqpminimizer_test_convergence(s, 1e-10, 1e-10);
      if(status == GSL_SUCCESS)
      {
         size_t j;
         if(verbose)
         {
            for(j=0;j<gsl_cqpminimizer_x(s)->size;j++)
            {

            }
         }
      }
      else{
         iter++;
      }
   } while (status == GSL_CONTINUE && iter<=max_iter);
   bool valid_solution = true;
   if(status != GSL_SUCCESS)
   {
      if(gsl_cqp_minimizer_test_infeasibility(s, 1e-10) != GSL_SUCCESS)
      {
         valid_solution = false;
      }
   }
   if(valid_solution)
   {
      gsl_vector_memcpy(x, gsl_cqpminimizer_x(s));
   }
   
}

void swapvals(int height, int width, PixelDoubleMatrix& one, PixelDoubleMatrix& two)
{
   for(int i = 0; i < height; i++)
   {
      for(int j=0; j < width; j++)
      {
         double tmp = one[i][j];
         one[i][j] = two[i][j];
         two[i][j] = tmp;
      }
   }
}
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
	pSrc->Convert(TMO_RGB); // This is format of Y as luminance
	pDst->Convert(TMO_RGB); // x, y as color information
   int imageHeight = pSrc->GetHeight();
   int imageWidth = pSrc->GetWidth();
	
												// of three doubles representing
												// three colour components

   //---------------------------------------Declarations of variables---------------------------------------------------------------------
   DisplaySize ds;
   DisplayFunc df;
   CPDfunction C;
   float res = 1024;
   float vd_screen = 2;
   float vd_size = 0.5f;
   float ref_white = 2;
   double cef = 1.f;
   float saturation = 1.f;
   displaySize(30.f, vd_size, ds);
   DisplayFuncInit(2.2f, 0.8, 200, 0.01, 60, df);
   CPDinit(C);
   double *imgSrcData = pSrc->GetData();
   double pixelR, pixelG, pixelB;
   float threshold = 0.0043;
   double adapt_scene = 1000;
   double stonits = pSrc->GetStonits();
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
         double tmp = rgb2luminance(pixelR, pixelG, pixelB);
         if(tmp < MIN_VALUE)
         {
            tmp = MIN_VALUE;
         }
         if(tmp > MAX_VALUE)
         {
            tmp = MAX_VALUE;
         }
         LogLuminancePixels[j][i] = log10(tmp);
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
            if(xi < 0 || xi >= C.x_cnt)
            {
               check_out = true;
               continue;
            }
            int gi = round_int((g + C.g_max)/C.delta);
            if(gi < 0 || gi >= C.g_cnt)
            {
               continue;
            }
            if(g > threshold && g < C.delta/2)
            {
               C.C_val[calcCval(xi,g_tp,i,C)] += 1.0;
            }
            else if(g < -1*threshold && g > -1*C.delta/2){
               C.C_val[calcCval(xi,g_tn,i,C)] += 1.0;
            }
            else{
               C.C_val[calcCval(xi,gi,i,C)] += 1.0;
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
         if(~grad)
         {
            C.C_val[calcCval(m,g_tp,i,C)] += 1.0;
            C.C_val[calcCval(m,g_tn,i,C)] += 1.0;
         }
         for(int k=0; k < C.g_cnt; k++)
         {
            if(k != g_t)
            {
               C.final_value += C.C_val[calcCval(m,k,i,C)];
            }
         }
      }
      swapvals(imageHeight,imageWidth,LogLuminancePixels, gausianOutVal);
   }
   fprintf(stderr,"final val %g\n",C.final_value);
   fprintf(stderr,"size %d\n",C.C_val.size());
   

   
	//------------------------------------------Computing tone curve-----------------------------------------------
   double display_drange = log10(calcDisplayFunc(1.f, df)/calcDisplayFunc(0.f, df));
   CSFvalues csf;
   ToneCurve tc;
   for(int f=0; f< C.freq_cnt; f++)
   {
      CSFvals tmp;
      tmp.delta = C.log_lum_scale[1] - C.log_lum_scale[0];
      tmp.v_size = C.x_cnt;
      //tmp.x_i = C.log_lum_scale;
      for(int i=0; i < C.x_cnt; i++)
      {
         tmp.y_i.push_back(cs_daly(C.f_scale[f],1,0,pow(10.,C.log_lum_scale[i])));
         tmp.x_i.push_back(C.log_lum_scale[i]);
      }
      csf.push_back(tmp);
   }
   int counter = 0;
   int max_ = (C.g_cnt - 1)/2;
   vector<int> used_v(C.x_cnt-1, 0);
   vector<int> unused(C.x_cnt-1);
   int min_max[2] = {C.x_cnt - 1, 0};

   for(int i=0; i < C.freq_cnt; i++)
   {
      for(int j=0; j < C.x_cnt; j++)
      {
         for(int k=max(0, j-max_); k < min(C.x_cnt-1, j + max_);k++)
         {
            if(i == j || C.C_val[calcCval(j, k-j+max_,i,C)] == 0)
            {
               continue;
            }
            counter++;
            int f = min(j,k);
            int t = max(j,k);
            for(int m=f; m < t; m++)
            {
               used_v[m] = 1;
            }
            min_max[0] = min(min_max[0], f);
            min_max[1] = max(min_max[1], t-1);
         }
      }
   }
   int white = 0;
   if(ref_white > 0)
   {
      counter++;
      float white_log = log10(ref_white);
      int i;
      for(i = C.x_cnt-1; i>= 0; i--)
      {
         if(C.log_lum_scale[i] <= white_log)
         {
            break;
         }
      }
      white = i;
      used_v[white] = 1;
      min_max[0] = min(min_max[0], white);
      min_max[1] = max(min_max[1], white);
   }

   int tmp = 0;
   int missing = 0;
   for(int m = 0; m < C.x_cnt-1; m++)
   {
      if(m < min_max[0] || m > min_max[1])
      {
         unused[m] = -1;
         continue;
      }
      if(!used_v[m])
      {
         if(m > 0 && !used_v[m-1])
         {
            unused[m] = -1;
            continue;
         }
         missing++;
      }
      unused[m] = tmp++;
   }
   int Eq_cnt = counter + missing;
   int Non_zero_var = tmp;
   gsl_matrix *A(gsl_matrix_calloc(Non_zero_var+1, Non_zero_var));
   gsl_matrix_set_identity(A);
   gsl_matrix_view lr = gsl_matrix_submatrix(A,Non_zero_var,0,1,Non_zero_var);
   gsl_matrix_set_all(&lr.matrix, -1);
   gsl_vector *D(gsl_vector_calloc(Non_zero_var+1));
   gsl_vector_set(D, Non_zero_var, -1*display_drange);

   gsl_matrix *M(gsl_matrix_calloc(Eq_cnt, Non_zero_var));
   gsl_vector *B(gsl_vector_alloc(Eq_cnt));
   gsl_vector *N(gsl_vector_alloc(Eq_cnt));

   size_t lum_background[Eq_cnt];
   size_t f_band[Eq_cnt];
   //for(int i=0; i < Eq_cnt; i++)
   //{
   //   lum_background[i] = 0;
   //}
   fprintf(stderr,"first counter %d\n",counter);
   counter = 0;
   for(int i=0; i < C.freq_cnt; i++)
   {
      double sens;
      if(adapt_scene != -1)
      {
         sens = cs_daly(C.f_scale[i],1,0,adapt_scene);
      }
      for(int j=0; j < C.x_cnt; j++)
      {
         for(int k=max(0, j - max_); k < min(C.x_cnt-1, j+max_);k++)
         {
            if(j==k || C.C_val[calcCval(j,k-j+max_,i,C)] == 0.0)
            {
               continue;
            }
            int f = min(j,k);
            int t = max(j,k);
            for(int m=f; m < t;m++)
            {
               if(unused[m] == -1)
               {
                  continue;
               }
               gsl_matrix_set(M, counter, unused[m], 1);
            }
            //TODO
            if(adapt_scene == -1)
            {
               sens = interpolation(C.log_lum_scale[f],csf[i]);
            }
            gsl_vector_set(B,counter, transducer((C.log_lum_scale[t] - C.log_lum_scale[f])*cef,sens));
            gsl_vector_set(N,counter, C.C_val[calcCval(j,k-j+max_,i,C)]);
            lum_background[counter] = j;
            f_band[counter] = i;
            counter++;
         }
      }
   }
   fprintf(stderr,"step\n");
   //TODO
   if(ref_white > 0)
   {
      for(int k=white; k < C.x_cnt-1; k++)
      {
         if(unused[k] == -1)
         {
            continue;
         }
         gsl_matrix_set(M, counter, unused[k], 1);
      }
      gsl_vector_set(B, counter, 0);
      gsl_vector_set(N, counter, C.final_value * 0.1);
      f_band[counter] = 0;
      lum_background[counter] = white;
      counter++;
   }
   for(int i = min_max[0]; i <= min_max[1]; i++)
   {
      if(!used_v[i])
      {
         int f = i;
         int t = i+1;
         while(!used_v[t])
         {
            t++;
         }
         assert(counter < Eq_cnt);
         for(int l=f; l < t; l++)
         {
            if(unused[l] == -1)
            {
               continue;
            }
            gsl_matrix_set(M, counter, unused[l],1);
         }
         double sens;
         if(adapt_scene == -1)
         {
            sens = interpolation(C.log_lum_scale[f],csf[C.freq_cnt-1]);
         }
         else
         {
            sens = cs_daly(C.f_scale[C.freq_cnt-1],1,0,adapt_scene);
         }
         gsl_vector_set(B,counter, transducer((C.log_lum_scale[t] - C.log_lum_scale[f])*cef,sens));
         gsl_vector_set(N,counter,C.final_value * 0.1);
         f_band[counter] = C.freq_cnt -1;
         lum_background[counter] = t;
         counter++;
         i = t;
      }
   }
   gsl_matrix *H(gsl_matrix_alloc(Non_zero_var, Non_zero_var));
   gsl_matrix *Na(gsl_matrix_alloc(Eq_cnt, Non_zero_var));
   gsl_matrix *Ak(gsl_matrix_alloc(Eq_cnt, Non_zero_var));
   gsl_vector *f(gsl_vector_alloc(Non_zero_var));
   gsl_vector *Ax(gsl_vector_alloc(Eq_cnt));
   gsl_vector *K(gsl_vector_alloc(Eq_cnt));
   gsl_vector *X(gsl_vector_alloc(Non_zero_var));
   gsl_vector *X_p(gsl_vector_alloc(Non_zero_var));
   for(int i=0; i < C.x_cnt; i++)
   {
      tc.y_i.push_back(0.0);
   }
   fprintf(stderr,"step2\n");
   //fprintf(stderr,"Eq: %d counter : %d\n",Eq_cnt,counter);
   //for(int i =0; i < Eq_cnt; i++)
   //{
   //   fprintf(stderr,"band: %d back: %d\n",f_band[i],lum_background[i]);
   //}
   gsl_vector_set_all(X, display_drange/Non_zero_var);
   int iterations = 200;
   for(int iter=0; iter < iterations; iter++)
   {
      calculateToneCurve(tc,unused,X,Non_zero_var,C.x_cnt,calcDisplayFunc(0,df),calcDisplayFunc(1,df));
      gsl_blas_dgemv(CblasNoTrans, 1, M, X, 0, Ax);
      //fprintf(stderr,"test1\n");
      for(int i=0; i < Eq_cnt; i++)
      {
         double tmp_var = gsl_vector_get(Ax, i);
         double t;
         if(fabs(tmp_var) < 0.0001)
         {
            t = 1.0;
         }
         else
         {
            t = tmp_var;
         }
         //for(int u=0; u < Eq_cnt; u++)
         //{
         //   fprintf(stderr,"%d : %d\n",u,lum_background[u]);
         //}
         //fprintf(stderr,"f %d x: %d\n",C.freq_cnt,C.x_cnt);
         //fprintf(stderr,"eq %d\n",Eq_cnt);
         //fprintf(stderr,"size %d index %d\n",tc.y_i.size(), lum_background[i]);
         //fprintf(stderr,"size2 %d index2 %d\n",csf.size(),f_band[i]);
         //fprintf(stderr,"xcnt %d eq %d\n",C.x_cnt, Eq_cnt); 
         int index = f_band[i];
         if(index >= csf.size() || index < 0)
         {
            index = csf.size()-1;
         }
         int index2 = lum_background[i];
         if(index2 >= tc.y_i.size() || index2 < 0)
         {
            index2 = tc.y_i.size()-1;
         }
         //fprintf(stderr,"lum: %d band: %d\n",lum_background[i],f_band[i]);
         double sens = interpolation(tc.y_i[index2],csf[index]);
         gsl_vector_set(K, i, transducer(tmp_var,sens)/t);
      }
      //fprintf(stderr,"test2\n");
      multiple(M,Ak,K);
      multiple(Ak, Na, N);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,Ak,Na,0,H);
      gsl_blas_dgemv(CblasTrans, -1, Na, B, 0, f);
      gsl_vector_memcpy(X_p, X);
      solver(H,f,A,D,X);

      double min_d = (C.log_lum_scale[1] - C.log_lum_scale[0])/10.0;
      bool c = true;
      for(int o=0; o < Non_zero_var; o++)
      {
         double del = fabs(gsl_vector_get(X, o) - gsl_vector_get(X_p, o));
         if(del > min_d)
         {
            c = false;
            break;
         }
      }
      if(c){
         break;
      }
   }
   calculateToneCurve(tc,unused,X,Non_zero_var,C.x_cnt,calcDisplayFunc(0,df),calcDisplayFunc(1,df));
   fprintf(stderr,"step3\n");
   double t_filter_a[] = { 1.000000000000000,  -2.748835809214676,   2.528231219142559,  -0.777638560238080 };
   double t_filter_b[] = { 0.000219606211225409,   0.000658818633676228,   0.000658818633676228,   0.000219606211225409 };
   //for(int k=0;k < tc.y_i.size();k++)
   //{
   //   fprintf(stderr,"%d : %g\n",k,tc.y_i[k]);
   //}

   CSFvals filterTC;
   for(int i=0; i < tc.y_i.size();i++)
   {
      filterTC.y_i.push_back(tc.y_i[i]);
   }
   for(int i=0; i < C.log_lum_scale.size();i++)
   {
      filterTC.x_i.push_back(C.log_lum_scale[i]);
   }
   
   /*for(int i=0; i < 4; i++)
   {
      for(int j=0; j < filterTC.y_i.size();j++)
      {
         filterTC.y_i[j] += t_filter_b[i] * tc.y_i[j];
         if(i > 0)
         {
            filterTC.y_i[j] -= t_filter_a[i] * tc.y_i[j];
         }
      }
   }*/
   for(int i=0; i < filterTC.y_i.size();i++)
   {
      if(filterTC.y_i[i] < log10(calcDisplayFunc(0,df)))
      {
         filterTC.y_i[i] = log10(calcDisplayFunc(0,df));
      }
      else if(filterTC.y_i[i] > log10(calcDisplayFunc(1,df)))
      {
         filterTC.y_i[i] = log10(calcDisplayFunc(1,df));
      }
      else
      {
         filterTC.y_i[i] = filterTC.y_i[i] ;
      }
      fprintf(stderr,"%d %g %g %g\n",i,filterTC.x_i[i], filterTC.y_i[i],calcInverseDisplayFunc((float)pow(10, filterTC.y_i[i]),df));
   }
   fprintf(stderr,"min: %g max: %g\n",log10(calcDisplayFunc(0.f,df)),log10(calcDisplayFunc(1.f,df)));








   //fprintf(stderr,"sizex %d sizetc %d\n",C.x_cnt, tc.y_i.size());
   CSFvals finalTC;
   finalTC.v_size = C.x_cnt;
   for(int i=0; i < filterTC.y_i.size();i++)
   {
      finalTC.y_i.push_back((float)pow(10, filterTC.y_i[i]));
   }
   for(int i=0; i < C.log_lum_scale.size();i++)
   {
      finalTC.x_i.push_back(C.log_lum_scale[i]);
   }
   finalTC.v_size = C.x_cnt;
   CSFvals cc;
   cc.v_size = C.x_cnt;
   for(int i=0; i < filterTC.y_i.size()-1; i++)
   {
      float contrast = std::max( (filterTC.y_i[i+1] - filterTC.y_i[i])/(C.log_lum_scale[i+1]-C.log_lum_scale[i]), 0.0 );
      float k1 = 1.48;
      float k2 = 0.82;
      cc.y_i.push_back(((1 + k1)*pow(contrast,k2))/(1 + k1*pow(contrast,k2)) * saturation );
      //fprintf(stderr,"CC %d: %g\n",i,cc.y_i[i]);
   }
   for(int i=0; i < finalTC.x_i.size();i++)
   {
      cc.x_i.push_back(finalTC.x_i[i]);
   }
   cc.y_i.push_back(1);
   //cc.v_size = C.x_cnt;
	
   fprintf(stderr,"test %g\n", display_drange/Non_zero_var);
   //pSrc->Convert(TMO_Yxy);
   double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
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
			//pY *= dParameter; // Parameters can be used like
							  // simple variables
         float tmp = rgb2luminance(pY,px,py);
         float l = minvalue(tmp);
         float lum = interpolation(log10(l), finalTC);
         float s = interpolation(log10(l), cc);
         //float s = 0.6;
         //fprintf(stderr,"%g   %f  %f   %f\n",tmp,l,lum,s);
			// and store results to the destination image
			*pDestinationData++ = calcInverseDisplayFunc(powf(minvalue(pY/l), s) * lum, df);
			*pDestinationData++ = calcInverseDisplayFunc(powf(minvalue(px/l), s) * lum, df);
			*pDestinationData++ = calcInverseDisplayFunc(powf(minvalue(py/l), s) * lum, df);
         //*pDestinationData++ = pow(minvalue(pY/l), saturation) * lum;
			//*pDestinationData++ = pow(minvalue(px/l), saturation) * lum;
			//*pDestinationData++ = pow(minvalue(py/l), saturation) * lum;
         //fprintf(stderr,"R: %g  G: %g  B: %g\n",minvalue(pY/l),minvalue(px/l),minvalue(py/l));
         //fprintf(stderr,"R: %g G: %g B: %g\n",calcInverseDisplayFunc(powf(minvalue(pY/l), s) * lum, df),calcInverseDisplayFunc(powf(minvalue(px/l), s) * lum, df),calcInverseDisplayFunc(powf(minvalue(py/l), s) * lum, df));
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
   //pDst->CorrectGamma(2.2);
	return 0;
}
