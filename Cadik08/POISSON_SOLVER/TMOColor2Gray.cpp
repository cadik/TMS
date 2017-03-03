//(c)Martin Cadik
//03/2007 - Prague

//POZOR - momentalne je to vypnuto (druhy return 0;)

//POZOR - neni dodelane:
//	double ZETA [psi_xmax][psi_ymax];
//	double PSI [psi_xmax][psi_ymax];
// -jsou "natvrdo" maximalni dimenze - pro mcdesk-test.hdr


//poisson solver based on CATAM,
//http://www.maths.cam.ac.uk/undergrad/tripos/catam/


/* --------------------------------------------------------------------------- *
 * TMOColor2Gray.cpp: implementation of the TMOColor2Gray class.                       *
 * --------------------------------------------------------------------------- */

const bool NEW_FORMULA=true;

#include "./TMOColor2Gray.h"

#undef CCATSLWIN

extern "C" {
	#include "catam.h"	
}

#include "RGB_Gamma.cpp"


/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOColor2Gray::TMOColor2Gray()
{
	SetName(L"Color2Gray");
	SetDescription(L"Color2Gray in gradient domain - version 0.1");

	normalize.SetName(L"normalize");
	normalize.SetDescription(L"Normalize output luminance (using 1/Lo_max) to be in <0, 1>");
	normalize.SetDefault(true);
	normalize=true;
	this->Register(normalize);

	gamma.SetName(L"gamma");
	gamma.SetDescription(L"Gamma Correction Factor [-]");
	gamma.SetDefault(2.2);
	gamma=2.2;
	gamma.SetRange(1.0e-3,1.0e+2);
	this->Register(gamma);

	s.SetName(L"s");
	s.SetDescription(L"overprojection step");
	s.SetDefault(1.);
	s=1.;
	s.SetRange(0.,2.);
	this->Register(s);
	
	eps.SetName(L"eps");
	eps.SetDescription(L"epsilon for gradient correction");
	eps.SetDefault(0.1);
	eps=0.1;
	eps.SetRange(0.,100.);
	this->Register(eps);
}

TMOColor2Gray::~TMOColor2Gray()
{
}


	const int psi_xmax=513;//2272;//400;
	const int psi_ymax=513;//2049;//513;
	double ZETA [psi_xmax][psi_ymax];
	double PSI [psi_xmax][psi_ymax];



double myPow(double num, double power){
	if(fabs(num)<EPS_8) return 0.;
	if(num>0) return (pow(num, power));
	else return(-pow(fabs(num), power));
}//myPow


double TMOColor2Gray::formula(double* data, long y1, long x1, long y2, long x2){
	const double wa=1./4.;
	const double wb=1./3.;

	long tmp_ind_1 = y1*xmax+x1;
	long tmp_ind_2 = y2*xmax+x2;

	double dL= data[3*tmp_ind_1]-data[3*tmp_ind_2];
	double da= data[3*tmp_ind_1+1]-data[3*tmp_ind_2+1];
	double db= data[3*tmp_ind_1+2]-data[3*tmp_ind_2+2];

	//return (dL);
	return (myPow( myPow(dL,3.) + myPow(wa*da,3.) + myPow(wb*db,3.), 1./3.));
	//return myPow(dL+w*da+w*db,1./3.);
	//return (dL+w*da+w*db);
	
}//formula



/* --------------------------------------------------------------------------- *
 * Gradient inconsistency correction                                           *
 * --------------------------------------------------------------------------- */
void TMOColor2Gray::inconsistencyCorrection(TMOImage& G_image, TMOImage& divG_image, double eps){
	long xmax=G_image.GetWidth(), ymax=G_image.GetHeight();
	long tmp_y, tmp_ind, i, j, iter=0;
	double *pG_image=G_image.GetData(), *pdivG_image=divG_image.GetData();
	double E=0, tmp_divG, maxE=0, avgE=0;

	do{
		avgE=maxE=0;
		for ( i = 0; i < ymax; i+=1 )
		{
			tmp_y = i*xmax;
			for ( j = 0; j < xmax; j+=1 )
			{
				tmp_ind=j+tmp_y;
				E=	pG_image[3*tmp_ind]		- //Gx
					pG_image[3*tmp_ind+1]	+ //Gy
					((j+1<xmax)?(pG_image[3*(tmp_ind+1)+1]):0) -	//Gy(i+1,j)
					((i+1<ymax)?(pG_image[3*(tmp_ind+xmax)]):0);	//Gx(i,j+1)

				if(fabs(E)>maxE) maxE=fabs(E);
				avgE+=fabs(E);

				E*=0.25*s;
				
				pG_image[3*tmp_ind]-=E;	//Gx
				pG_image[3*tmp_ind+1]+=E;	//Gy
				if(j+1<xmax) pG_image[3*(tmp_ind+1)+1]-=E;		//Gy(i+1,j)
				if(i+1<ymax) pG_image[3*(tmp_ind+xmax)]+=E;	//Gx(i,j+1)
				
				tmp_divG=(((j-1)<0)?(pG_image[3*tmp_ind]):(pG_image[3*tmp_ind]-pG_image[3*(tmp_ind-1)]));
				tmp_divG=tmp_divG+(((i-1)<0)?(pG_image[3*tmp_ind+1]):(pG_image[3*tmp_ind+1]-pG_image[3*(tmp_ind-xmax)+1]));//*/				
				pdivG_image[3*tmp_ind]=tmp_divG;
				pdivG_image[3*tmp_ind+1]=tmp_divG;
				pdivG_image[3*tmp_ind+2]=tmp_divG;
			}
		}
		iter++;
		if(!(iter%100)){
			printf("Iteration: %d,\t maxE=%8.2g, \tavgE=%8.2g \n", iter, maxE, avgE);
		}		
	}while(maxE>eps);
}//inconsistencyCorrection


/* --------------------------------------------------------------------------- *
 * Gradient field integration                                                  *
 * --------------------------------------------------------------------------- */
void TMOColor2Gray::GFintegration(TMOImage& G_image, TMOImage& Dst_image)
{
	long xmax=Dst_image.GetWidth(), ymax=Dst_image.GetHeight();
	long tmp_y, tmp_ind, i, j;
	double *pG_image=G_image.GetData(), *pDst_image=Dst_image.GetData();

	pDst_image[0]=pDst_image[1]=pDst_image[2]=0.;
	for ( i = 0; i < ymax; i++ )
	{
		tmp_y = i*xmax;
		if(i>0) pDst_image[3*tmp_y]=pDst_image[3*tmp_y+1]=pDst_image[3*tmp_y+2]=
			(pDst_image[3*(tmp_y-xmax)] + pG_image[3*(tmp_y-xmax)+1]);
			//neboli: OUTPUT_BW[0][y] = OUTPUT_BW[0][y-1] + Grad_Y[0][y-1];

		for ( j = 1; j < xmax; j++ )
		{
			tmp_ind=j+tmp_y;
			pDst_image[3*tmp_ind]=pDst_image[3*tmp_ind+1]=pDst_image[3*tmp_ind+2]=
				(pDst_image[3*(tmp_ind-1)] + pG_image[3*(tmp_ind-1)]);
			//neboli: OUTPUT_BW[x][y] = OUTPUT_BW[x-1][y] + Grad_X[x-1][y]; 
		}
	}
}//GFintegration



/* --------------------------------------------------------------------------- *
 * Calibration of the output values                                            *
 * --------------------------------------------------------------------------- */
void TMOColor2Gray::calibrate(TMOImage& src_image, TMOImage& dst_image){
	long i, j, tmp_y, xmax=src_image.GetWidth(), ymax=src_image.GetHeight();
	double	*pSrc_image=src_image.GetData(),
			*pDst_image=dst_image.GetData();
	double  SUM_L_new=0, SUM_L_old=0, SUM_L2_new=0, SUM_L_newL_old=0;
	double	A=0, B=0;

	//assert: src_image je v Yxy
	//assert: dst_image je v RGB a to v stupnich sedi

	for ( i = 0; i < ymax ; i++ )
	{
		tmp_y = i*xmax;
		for ( j = 0; j < xmax ; j++ )
		{
			SUM_L_new+=pDst_image[3*(tmp_y+j)];
			SUM_L_old+=pSrc_image[3*(tmp_y+j)];
			SUM_L2_new+=pDst_image[3*(tmp_y+j)]*pDst_image[3*(tmp_y+j)];
			SUM_L_newL_old+=pDst_image[3*(tmp_y+j)]*pSrc_image[3*(tmp_y+j)];
			
		}
	}

	if(SUM_L2_new - SUM_L_new * SUM_L_new!=0)
		B=(SUM_L_newL_old - SUM_L_new * SUM_L_old)/(SUM_L2_new - SUM_L_new * SUM_L_new);
	A=(SUM_L_old - B*SUM_L_new)/(xmax*ymax); 
	printf("Normalization: A+B*L == %g+%g*L\n", A, B);

	for ( i = 0; i < ymax ; i++ )
	{
		tmp_y = i*xmax;
		for ( j = 0; j < xmax ; j++ )
		{
			pDst_image[3*(tmp_y+j)]=A+B*pDst_image[3*(tmp_y+j)];
			pDst_image[3*(tmp_y+j)+1]=A+B*pDst_image[3*(tmp_y+j)+1];
			pDst_image[3*(tmp_y+j)+2]=A+B*pDst_image[3*(tmp_y+j)+2];
		}
	}


}//calibrate




/* --------------------------------------------------------------------------- *
 * An implementation of Poisson solver                                         *
 * --------------------------------------------------------------------------- */
TMOColor2Gray::Transform()
{
	int i=0, j=0;
	double L_max=0.,  // max luminance 
		   L_min=0.,  // min luminance
		   L_temp=0.;   // 
	double dStonits=pSrc->GetStonits();
	pSrc->Convert(TMO_Yxy);
	pSrc->GetMinMaxAvg(&L_min, &L_max, &L_temp);
	fprintf(stdout, "Min luminance: %g[cd/m^2]\nMax luminance: %g[cd/m^2]\nAvg luminance: %g[cd/m^2]\n", L_min * dStonits, L_max * dStonits, L_temp * dStonits);
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB, true);
	xmax=pSrc->GetWidth();
	ymax=pSrc->GetHeight();
	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();
	int max,k;
	int tmp_y,tmp_ind;


///////////////////////////////////////////////////////////////////////////////////////
// Vypocet gradientu

	int max_square_pow2=MAX(xmax, ymax);
	double log2=log((double)max_square_pow2) / log(2.0);
	if(max_square_pow2 > pow(2.0, floor(log2))) 
		max_square_pow2=pow(2.0, ceil(log2));
	max=max_square_pow2*max_square_pow2;
	int jmax = max_square_pow2;
	jmax++;

	double **H=alloc_image1(max_square_pow2, max_square_pow2);
	multiply_image1(H, max_square_pow2, max_square_pow2, 0.0); 
	vect2D *nablaH=new vect2D[max];
	double tmp_divG;
	TMOImage G_image;
	TMOImage divG_image;

	if(!NEW_FORMULA){
	/////////////////////////////////////////////////////////////
	//////vypocet gradientu z luminance
	//// vypocte hodnoty luminance pro vsechny pixely a zaroven log luminance (H)
	pSrc->Convert(TMO_Yxy);
	for (i = 0; i < ymax; i++)
	{
		for (int j = 0; j < xmax; j++)
		{
			double L_w = *pSourceData++;
			double x_w = *pSourceData++;
			double y_w = *pSourceData++;

			H[i][j]=TAKE_LOG(L_w);
			//H[i][j]=L_w;
	
			*pDestinationData++ = L_w ;  
			*pDestinationData++ = x_w; //colors remain unchanged
			*pDestinationData++ = y_w; //colors remain unchanged
		}
	}	
	
	// vypocte hodnoty gradientu H (nabla H)
	// a odhadne parametr alpha jako 0.1 * Avg ||\Nabla H_k(x, y)||
	double avg_gradient=0;
	for (i = 0; i < ymax; i++)
	{
		tmp_y = i*xmax;
		for (j = 0; j < xmax; j++)
		{
			tmp_ind=j+tmp_y;
			nablaH[tmp_ind].x=(((j+1)==xmax)?0.0:(H[i][j+1]-H[i][j]));
			nablaH[tmp_ind].y=(((i+1)==ymax)?0.0:(H[i+1][j]-H[i][j]));
			avg_gradient+=sqrt(nablaH[tmp_ind].x*nablaH[tmp_ind].x + nablaH[tmp_ind].y*nablaH[tmp_ind].y);
		}
	}
	avg_gradient/=(xmax*ymax);
	free_image1(H, max_square_pow2, max_square_pow2);
	
	// vypocteme hodnoty divG
	//a taky si je ulozime G a divG -> dump:
	G_image.New(xmax, ymax); //Gx,Gy,DivG
	divG_image.New(xmax, ymax);
	double *pG_image = G_image.GetData();
	double *pdivG_image = divG_image.GetData();
	int image_tmp_y=0, image_tmp_ind=0;
	char filename[500];
	strcpy(filename, pSrc->GetFilename());

	for ( i = 0; i < ymax; i++ )
	{
		tmp_y = i*xmax;
		image_tmp_y = i*jmax;
		for ( j = 0; j < xmax; j++ )
		{
			tmp_ind=j+tmp_y;
			image_tmp_ind=j+image_tmp_y;
			tmp_divG=(((j-1)<0)?(nablaH[tmp_ind].x):(nablaH[tmp_ind].x-nablaH[tmp_ind-1].x));
			tmp_divG=tmp_divG+(((i-1)<0)?(nablaH[tmp_ind].y):(nablaH[tmp_ind].y-nablaH[(i-1)*xmax+j].y));//*/

			pG_image[3*tmp_ind]=fabs(nablaH[tmp_ind].x); //Gx
			pG_image[3*tmp_ind+1]=fabs(nablaH[tmp_ind].y); //Gy
			pG_image[3*tmp_ind+2]=0; //

			pdivG_image[3*tmp_ind]=tmp_divG;
			pdivG_image[3*tmp_ind+1]=tmp_divG;
			pdivG_image[3*tmp_ind+2]=tmp_divG;
		}
	}

	G_image.SetFilename(filename);
	G_image.SaveWithSuffix("G");
	//G_image.SaveWithSuffix("G", TMO_RAW);	

	divG_image.SetFilename(filename);
	divG_image.SaveWithSuffix("divG");
	//divG_image.SaveWithSuffix("divG", TMO_RAW);	
	}
	else
	{
	/////////////////////////////////////////////////////////////
	//////vypocet gradientu z Lab
	//prevod do Lab
	double X, Y, Z;
	for (i = 0; i < ymax; i++)
	{
		tmp_y = i*xmax;
		for (j = 0; j < xmax; j++)
		{
			tmp_ind = tmp_y+j;

			RGB709_XYZ(Cdisplay_01_Clinear(pSourceData[3*tmp_ind]), Cdisplay_01_Clinear(pSourceData[3*tmp_ind+1]), Cdisplay_01_Clinear(pSourceData[3*tmp_ind+2]), &X, &Y, &Z);
			XYZ_Lab_(X, Y, Z, Xn, Yn, Zn, &pSourceData[3*tmp_ind], &pSourceData[3*tmp_ind+1], &pSourceData[3*tmp_ind+2]);
		}
	}

	//// vypocte hodnoty gradientu H (nabla H)
	double avg_gradient=0;
	for (i = 0; i < ymax; i++)
	{
		tmp_y = i*xmax;
		for (j = 0; j < xmax; j++)
		{
			tmp_ind=j+tmp_y;
			nablaH[tmp_ind].x=(((j+1)==xmax)?0.0:formula(pSourceData, i,j+1,i,j));
				//(H[i][j+1]-H[i][j]));
			nablaH[tmp_ind].y=(((i+1)==ymax)?0.0:formula(pSourceData, i+1,j,i,j));
				//(H[i+1][j]-H[i][j]));
			avg_gradient+=sqrt(nablaH[tmp_ind].x*nablaH[tmp_ind].x + nablaH[tmp_ind].y*nablaH[tmp_ind].y);
		}
	}
	avg_gradient/=(xmax*ymax);
	free_image1(H, max_square_pow2, max_square_pow2);


	// vypocteme hodnoty divG
	//a taky si je ulozime G a divG -> dump:
	G_image.New(xmax, ymax); //Gx,Gy,DivG
	divG_image.New(xmax, ymax);
	double *pG_image = G_image.GetData();
	double *pdivG_image = divG_image.GetData();
	int image_tmp_y=0, image_tmp_ind=0;
	char filename[500];
	strcpy(filename, pSrc->GetFilename());

	for ( i = 0; i < ymax; i++ )
	{
		tmp_y = i*xmax;
		image_tmp_y = i*jmax;
		for ( j = 0; j < xmax; j++ )
		{
			tmp_ind=j+tmp_y;
			image_tmp_ind=j+image_tmp_y;
			tmp_divG=(((j-1)<0)?(nablaH[tmp_ind].x):(nablaH[tmp_ind].x-nablaH[tmp_ind-1].x));
			tmp_divG=tmp_divG+(((i-1)<0)?(nablaH[tmp_ind].y):(nablaH[tmp_ind].y-nablaH[(i-1)*xmax+j].y));//*/

			pG_image[3*tmp_ind]=nablaH[tmp_ind].x; //Gx
			pG_image[3*tmp_ind+1]=nablaH[tmp_ind].y; //Gy
			pG_image[3*tmp_ind+2]=0; //

			pdivG_image[3*tmp_ind]=tmp_divG;
			pdivG_image[3*tmp_ind+1]=tmp_divG;
			pdivG_image[3*tmp_ind+2]=tmp_divG;
		}
	}

	G_image.SetFilename(filename);
	G_image.SaveWithSuffix("G");
	G_image.SaveWithSuffix("G", TMO_EXR);	

	divG_image.SetFilename(filename);
	divG_image.SaveWithSuffix("divG");
	//divG_image.SaveWithSuffix("divG", TMO_RAW);	

	inconsistencyCorrection(G_image, divG_image, eps);
	G_image.SetFilename(filename);
	G_image.SaveWithSuffix("G_corrected");
	divG_image.SetFilename(filename);
	divG_image.SaveWithSuffix("divG_corrected");

	GFintegration(G_image, *pDst);

	double *pDst_image=pDst->GetData();
	for ( i = 0; i < pDst->GetHeight() ; i++ )
	{
		tmp_y = i*(pDst->GetWidth());
		for ( j = 0; j < pDst->GetWidth() ; j++ )
		{
			pDst_image[3*(tmp_y+j)]+=100;
			pDst_image[3*(tmp_y+j)+1]+=100;
			pDst_image[3*(tmp_y+j)+2]+=100;
		}
	}

	pDst->Convert(TMO_Yxy);
	double L_temp2;
	pDst->GetMinMaxAvg(&L_min, &L_max, &L_temp2);
	//fprintf(stdout, "Min luminance: %g[cd/m^2]\nMax luminance: %g[cd/m^2]\nAvg luminance: %g[cd/m^2]\n", L_min * dStonits, L_max * dStonits, L_temp * dStonits);


	pDst->Convert(TMO_RGB);

	calibrate(*pSrc, *pDst);

	//for ( i = 0; i < pDst->GetHeight() ; i++ )
	//{
	//	tmp_y = i*(pDst->GetWidth());
	//	for ( j = 0; j < pDst->GetWidth() ; j++ )
	//	{
	//		pDst_image[3*(tmp_y+j)]*=1/L_max;
	//		pDst_image[3*(tmp_y+j)+1]*=1/L_max;
	//		pDst_image[3*(tmp_y+j)+2]*=1/L_max;
	//	}
	//}


	pDst->Convert(TMO_Yxy);
	pDst->GetMinMaxAvg(&L_min, &L_max, &L_temp);
	fprintf(stdout, "Min luminance: %g[cd/m^2]\nMax luminance: %g[cd/m^2]\nAvg luminance: %g[cd/m^2]\n", L_min * dStonits, L_max * dStonits, L_temp * dStonits);

	//pDst->Convert(TMO_RGB);
	//pDst->CorrectGamma(gamma);

	pSrc->Convert(TMO_RGB);
	pSrc->SaveWithSuffix("INPUT-test");
	return 0;

	}





///////////////////////////////////////////////////////////////////////////////////////
// POISSON SOLVER
//	
// Poisson - priprava poli 
	long nx=xmax;
	long ny=max_square_pow2;
	double dlx=1;//0.001953125;//0.1;
	double dly=1;//0.001953125;//0.125;
	long psimax=xmax*max_square_pow2*(2*10);

	//double* psi= new double[psimax];
	//double* zeta= new double[psimax];
	//double* ZETA = new double[psi_xmax*psi_ymax];
	//double* PSI = new double[psi_xmax*psi_ymax];
	//double** ZETA=alloc_image1(psi_xmax, psi_ymax);
	//double** PSI=alloc_image1(psi_xmax, psi_ymax);


//	for (i=0;i<psimax;i++)
//	{
//		psi[i]=0.;
//		zeta[i]=0.;		
//	}

/*debug*/  for(i=0; i<psi_xmax;i++)
/*debug*/		for(j=0; j<psi_ymax; j++)
/*debug*/		{	ZETA[i][j]=0;
/*debug*/			PSI[i][j]=0;}

/*debug*/	TMOImage zetaI;
/*debug*/	zetaI.New(nx, ny, TMO_Y);
	double tmp_div=0;
	double *pdivG=divG_image.GetData();
	double *pzeta=zetaI.GetData();
	tmp_ind=0;

	int shift_i=10, shift_j=10;
	//int shift_i=0, shift_j=0;


	for (i = 0; i < xmax+2*shift_i; i++)
	{
		for (j = 0; j < ymax+2*shift_j; j++){
			ZETA[i][j]=0;
		}
	}	


	for (i = 0; i < xmax; i++)
	{
		for (j = 0; j < ymax; j++){
			tmp_div=*pdivG++; pdivG+=2;
//			zeta[tmp_ind++]=tmp_div;
/*debug*/	*pzeta++=tmp_div;
			ZETA[shift_i+i][shift_j+j]=divG_image.GetPixel(i, j)[0];
		}
	}	
/*debug*/	zetaI.SaveAs("zeta.tiff");  



//	inconsistencyCorrection(ZETA, psi_xmax, psi_ymax-1);
	pzeta=zetaI.GetData();
	for (i = 0; i < xmax; i++)
	{
		for (j = 0; j < ymax; j++){
			tmp_div=*pdivG++; pdivG+=2;
			*pzeta++=ZETA[i][j];
		}
	}	
/*debug*/	zetaI.SaveAs("zeta_inconsistencyCorrected.tiff");  



/////////////////////////////////////////////////////////////////////////////////////
// reseni Poissonovy rovnice a prekopirovani do vysledneho obr.
	printf("Solving Poisson equation....\n");
	//PoissonCL(psi, zeta, nx, ny, dlx, dly);
	//PoissonCL((two_d_array)PSI, (two_d_array)ZETA, nx, ny, dlx, dly);
	PoissonCL((two_d_array)PSI, (two_d_array)ZETA, psi_xmax, psi_ymax-1, 1, 1);

/*debug*/	TMOImage zetI;
/*debug*/	zetI.New(psi_xmax, psi_ymax, TMO_Y);
/*debug*/	for(i=0; i<psi_xmax; i++)
/*debug*/		for(j=0; j<psi_ymax; j++)
/*debug*/			*zetI.GetPixel(i,j)=(ZETA[i][j]);
/*debug*/	zetI.SaveAs("zet.tiff");
	
/*debug*/	TMOImage psiI;
/*debug*/	psiI.New(psi_xmax, psi_ymax, TMO_Y);
/*debug*/	for(i=0; i<psi_xmax; i++)
/*debug*/		for(j=0; j<psi_ymax; j++)
/*debug*/			*psiI.GetPixel(i,j)=(PSI[i][j]);
/*debug*/	psiI.SaveAs("psi.tiff");


////////////////////////////////////////////////////////////////////////////////////////
// presun z vysledneho pole, desaturace
	double sumR=0, sumG=0, sumB=0;
	double L_out, min_lum=HUGE_VAL, max_lum=-HUGE_VAL, avg_lum=0;

	pDestinationData=pDst->GetData();
	for (i = 0; i < ymax; i++)
	{
		tmp_y = i*xmax;
		for (j = 0; j < xmax; j++)
		{
			tmp_ind = tmp_y+j;		
			if(NEW_FORMULA){
				max_lum=MAX(max_lum, PSI[j+shift_i][i+shift_j]);
				min_lum=MIN(min_lum, PSI[j+shift_i][i+shift_j]);
				avg_lum+=L_out;
			}
			else
			{
				max_lum=MAX(max_lum, exp(PSI[j+shift_i][i+shift_j]));
				min_lum=MIN(min_lum, exp(PSI[j+shift_i][i+shift_j]));
				avg_lum+=L_out;
			}
		}
	}
	printf("Poisson solver output values range: [%f, %f]", min_lum, max_lum);
	avg_lum/=(xmax*ymax);
	printf("\nAvg output: %f\n", avg_lum);

	double R,G,B;
	for (i = 0; i < ymax; i++)
	{
		tmp_y = i*xmax;
		for (j = 0; j < xmax; j++)
		{			
			tmp_ind = tmp_y+j;
			if(NEW_FORMULA)
			{
			////truncation [0,100]
			//if(PSI[j][i]<0) PSI[j][i]=0;
			//else if(PSI[j][i]>100) PSI[j][i]=100;
			//shifting
			PSI[j+shift_i][i+shift_j]+=fabs(min_lum);
			CIE_Y_from_L(PSI[j+shift_i][i+shift_j], &L_out);
			RGBlinear_RGBdisplay(L_out, L_out, L_out, &R, &G, &B);
			pDestinationData[3*tmp_ind]=R; 
			pDestinationData[3*tmp_ind+1]=G;
			pDestinationData[3*tmp_ind+2]=B;
			//pDestinationData[3*tmp_ind]=(L_out+fabs(min_lum))/100.;
			//pDestinationData[3*tmp_ind+1]=(L_out+fabs(min_lum))/100.;
			//pDestinationData[3*tmp_ind+2]=(L_out+fabs(min_lum))/100.;			
			
			pDestinationData[3*tmp_ind]=L_out;
			pDestinationData[3*tmp_ind+1]=L_out;
			pDestinationData[3*tmp_ind+2]=L_out;			
			}	
			else{
				L_out=PSI[j+shift_i][i+shift_j];
				pDestinationData[3*tmp_ind]=exp(L_out);
				pDestinationData[3*tmp_ind+1]=exp(L_out);
				pDestinationData[3*tmp_ind+2]=exp(L_out);			
			}
		}
	}
	

/////////////////////////////////////////////////////////////////////////////////////////////
//// normalizace, gamma korekce
	//normalization
	if(normalize){
		max_lum=fabs(max_lum);
		printf("Normalization of output luminance by factor 1/%g to be in <0, 1>\n", max_lum);
		for (i = 0; i < ymax; i++)
		{
			tmp_y = i*xmax;
			for (j = 0; j < xmax; j++)
			{
				tmp_ind = 3*(tmp_y+j);
				//normalization
				pDestinationData[tmp_ind]/=max_lum;
				pDestinationData[tmp_ind+1]/=max_lum;
				pDestinationData[tmp_ind+2]/=max_lum;
			}
		}
	}//normalize

	// gamma correction
	pDst->Convert(TMO_RGB);
	if(NEW_FORMULA) pDst->CorrectGamma(gamma);
	pDst->ProgressBar(1, 1);

	return 0;
}//Transform
