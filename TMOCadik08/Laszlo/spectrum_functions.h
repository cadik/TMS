//spectrum functions: spectrum -> XYZ
//1nm in [360,830] interval, CIE 1931
//D65 given in [300,830]
//


char jel, jelsor[800];

double _read_data(FILE * file)
{

double szam;

long k, volt_szamjegy;

for(k=0;k<800;k++) jelsor[k] = 0;

k = volt_szamjegy = 0;
	
while(1)
	{
		jel = (unsigned char) fgetc(file);
		jelsor[k] = jel;
		k++;
		if('0'<=jel && jel<='9') volt_szamjegy++;
		if(   volt_szamjegy > 0 &&
			 !(jel == '.' || (jel >= '0' && jel <= '9') || jel == '-' || jel == '+' )
		  ) 
		break;
	}

szam = atof(jelsor);

return(szam);

}

BYTE READ_METAMER(long * szamlal)
{
		FILE * in_metamer;
		long k;
		double X, Y, Z;

		if((in_metamer = fopen("metamer_in.txt","rt"))==0)
			return(0);
		
	*szamlal = 0;

ciklus:
	k = (long) _read_data(in_metamer);
	printf("\n\t --- tipus = %ld",k);

	metamer[*szamlal].tipus = k;

	if(k<0 || k>2) goto close;

		X = _read_data(in_metamer);
		printf("\n\t elso     %g",X);
		Y = _read_data(in_metamer);
		printf("\n\t masodik  %g",Y);
		Z = _read_data(in_metamer);
		printf("\n\t harmadik %g",Z);
		metamer[*szamlal].X = X;
		metamer[*szamlal].Y = Y;
		metamer[*szamlal].Z = Z;
		(*szamlal)++;
		goto ciklus;

close:
		fclose(in_metamer);
return(1);
}

void READ_SPECTRUM_DATAE(void)
{

	FILE * spectrum_in;
	short i;
	double s,t;
	//char c;

	spectrum_in   =  fopen("Cadik08/Laszlo/Xyz31_1.txt",  "rt");

	for(i=360;i<=830;i++)
	{
		//fscanf(spectrum_in,"%d,%f,%f,%f",&vawelength,&x,&y,&z);

		_read_data(spectrum_in);
		CIE1931_X[i] = _read_data(spectrum_in);
		CIE1931_Y[i] = _read_data(spectrum_in);
		CIE1931_Z[i] = _read_data(spectrum_in);

		//ezek 100.-szorosai a hatarszinek XYZ koordinatai

		t = CIE1931_X[i] + CIE1931_Y[i] + CIE1931_Z[i];
		t = 1./t;

		CIE1931_x[i] = t * CIE1931_X[i];
		CIE1931_y[i] = t * CIE1931_Y[i];
		CIE1931_z[i] = t * CIE1931_Z[i];

			
		t = (CIE1931_y[i]-D65_y) / (CIE1931_x[i]-D65_x);
		fi_D65[i]=RAD_GRAD * atan2( (CIE1931_y[i]-D65_y), (CIE1931_x[i]-D65_x) );

		fi_C_[i]=RAD_GRAD * atan2(CIE1931_y[i]-C_y, CIE1931_x[i]-C_x);

	}

	fclose(spectrum_in);

	spectrum_in = fopen("Cadik08/Laszlo/LIGHTD65.TXT","rt");

	for(i=300;i<=830;i++)
	{
		_read_data(spectrum_in);
		light_D65[i] = _read_data(spectrum_in);
	}
	
	fclose(spectrum_in);

	for(i=360,s=0;i<=830;i++)
		s += light_D65[i] * CIE1931_Y[i];

		t = 100. / s;

	for(i=360;i<=830;i++)
	light_D65_normalized[i] = light_D65[i] * t;

	//for(i=360,s=0;i<=830;i++)
	//	s += light_D65_normalized[i] * CIE1931_y[i];
	//ennek 100 az Y erteke konstans 1 reflektancia mellett

	spectrum_in = fopen("Cadik08/Laszlo/light_A.txt","rt");

	for(i=300;i<=830;i++)
	{
		_read_data(spectrum_in);
		light_A[i] = _read_data(spectrum_in);        
	}

	for(i=360,s=0;i<=830;i++)
		s += light_A[i] * CIE1931_Y[i];

		t = 100. / s;

	for(i=360;i<=830;i++)
	light_A_normalized[i] = light_A[i] * t;
	
	fclose(spectrum_in);
}

void spectrum_1nm_D65_XYZ(double *spectrum,
			double *X,  double *Y,  double *Z,
			double *Xn, double *Yn, double *Zn)
{
//=========================================
//CIE 1931 for 2 grad view field

	long i;

	*X = *Y = *Z = *Xn = *Yn = *Zn = 0;

	for(i=360;i<831;i++)
	{
		*Xn += light_D65[i]*CIE1931_X[i];
		*Yn += light_D65[i]*CIE1931_Y[i];
		*Zn += light_D65[i]*CIE1931_Z[i];

		*X  += light_D65[i]*CIE1931_X[i]*spectrum[i];
		*Y  += light_D65[i]*CIE1931_Y[i]*spectrum[i];
		*Z  += light_D65[i]*CIE1931_Z[i]*spectrum[i];
	}

}

void spectrum_1nm_A_light_XYZ(double *spectrum,
			double *X,  double *Y,  double *Z,
			double *Xn, double *Yn, double *Zn)
{
//=========================================
//CIE 1931 for 2 grad view field

	long i;

	*X = *Y = *Z = *Xn = *Yn = *Zn = 0;

	for(i=360;i<831;i++)
	{
		*Xn += light_A[i]*CIE1931_X[i];
		*Yn += light_A[i]*CIE1931_Y[i];
		*Zn += light_A[i]*CIE1931_Z[i];

		*X  += light_A[i]*CIE1931_X[i]*spectrum[i];
		*Y  += light_A[i]*CIE1931_Y[i]*spectrum[i];
		*Z  += light_A[i]*CIE1931_Z[i]*spectrum[i];
	}

}

double PLANCK(double T, double lambda)
{
// lambda nanometerben (celszeruen 360...830 intervallumban,
// T a fekete test hofoka Kelvinben,
// ld.: Stiles - Wyszecky,2nd ed., 1982, pp.12-13.

	double L,M,s;

        L = 1E-9 * lambda;
        s = c2 / L / T;
        M = 1e30 * c1 * pow(1E-3*lambda, -5.)/(exp(s) - 1.);
        // [M] = [Wm^(-3)
		
        return(M);
}

void spectrum_1nm_Planck_T_XYZ(double T, double *spectrum,
					double *X,  double *Y,  double *Z,
					double *Xn, double *Yn, double *Zn)
{
//========================================================

	double s;
	long i;

	*X = *Y = *Z = *Xn = *Yn = *Zn = 0;

	for(i=360;i<831;i++)
	{
		s = PLANCK(T,(double)i);

		*Xn += s * CIE1931_X[i];
		*Yn += s * CIE1931_Y[i];
		*Zn += s * CIE1931_Z[i];

		*X  += s * CIE1931_X[i] * spectrum[i];
		*Y  += s * CIE1931_Y[i] * spectrum[i];
		*Z  += s * CIE1931_Z[i] * spectrum[i];
	}

}

void Planck_series(void)
{
long i;
double s, T, X, Y, Z, Xn, Yn, Zn;

	for(i=1;i<831;i++)
		spectrum[i] = 1;

	for(i=0;i<20;i++)
	{
		T = 1500. + i * 500.;
		spectrum_1nm_Planck_T_XYZ(T, spectrum,
					&X, &Y, &Z,	&Xn, &Yn, &Zn);
	
		s = 100. / Y;
		X *= s;
		Z *= s;
		Y = 100.;

		printf("\n\n\tPlanck T = %6.1lf K,   X = %5.2lf   Y = %5.2lf   Z = %5.2lf",
					T, X, Y, Z);
	}

}




double max_tel1(long i, double V)
{

	double Y, Y_limit, T_max;
	T_max = 0.;

	Y = 0.01 * V * V;
	Y_limit = MAX_TEL[i].Y;

	if(Y >= Y_limit) T_max = (100. - Y) / (100. - Y_limit) * MAX_TEL[i].T;
	else             T_max =         Y  /         Y_limit  * MAX_TEL[i].T;

	return(T_max);
}



int READ_COLOR2GRAY_DATAE(void)
{
	FILE *inp, *tab;
	long i, j, k, i1, j1, tel;
	double g, s, _A, _T, _V, _V45, _V65, _V85, _T45, _T65, _T85;

        inp =  fopen("Cadik08/Laszlo/COLOR_GRAY_UNIF.txt","rt");
        tab =  fopen("Cadik08/Laszlo/GOLOR2GRAY_datae.txt","wt");


		for(i = 1; i <= 16; i++)
		{
			parok_rel_gray[i] = _read_data(inp);
			parok[i].A1 = _read_data(inp);
		    parok[i].T1 = _read_data(inp);
		    parok[i].V1 = _read_data(inp);

			parok[i].A2 = _read_data(inp);
		    parok[i].T2 = _read_data(inp);
		    parok[i].V2 = _read_data(inp);
			g = _read_data(inp);

			if(fabs(g + parok_rel_gray[i]) > 1E-8) return(1);

			fprintf(tab,"\n\n\ti=%2ld\trel szurke =  %5.1lf\tA1=%5.1lf\tT1=%5.1lf\tV1=%5.1lf\tA2=%5.1lf\tV2=%5.1lf\tT2=%5.1lf",
					     i,parok_rel_gray[i],parok[i].A1,parok[i].T1,parok[i].V1,parok[i].A2,parok[i].V2,parok[i].T1);
				

		}

		//--------------------------------------------------------------------------------------------------

		for(i=1; i <=7; i++)
			for(j=1; j <=7; j++)
			{
				_read_data(inp);
				_read_data(inp);
				HUE_PAIRS[i][j] = s = _read_data(inp);
				if(s < 0) sig[i][j] = -1; else sig[i][j] = +1;
				if(i == j) HUE_PAIRS[i][j] = 0.;
			}
		
		for(i=1; i <=7; i++)
			for(j=1; j < i; j++)
			{
				s = fabs(HUE_PAIRS[i][j])+fabs(HUE_PAIRS[j][i]);
				s /= 2.;
				HUE_PAIRS[i][j] = s * sig[i][j];
				HUE_PAIRS[i][j] = - HUE_PAIRS[j][i];

				//fuggveny hue parokra a T_rel relativ telitettsegek 0...5 es azon t u l i  ertekeire
				//x = T_rel < 1:   0.5 * x,  else   sqrt(x) - 0.5;  ertek es derivaltfolytonos ************************
			}

		for(i=1; i <=7; i++)
		{
			fprintf(tab,"\n");
			//for(j=1; j <=7; j++)
			//	fprintf(tab,"\n\t\tA1=%2ld\t\tA2=%2ld\t\trel szurke =  %5.2lf",
			//		10*i,10*j,HUE_PAIRS[i][j]);
			
				fprintf(tab,"\n\t\tA = %2ld\t%5.2lf\t%5.2lf\t%5.2lf\t%5.2lf\t%5.2lf\t%5.2lf\t%5.2lf",
					10*i,HUE_PAIRS[i][1],HUE_PAIRS[i][2],HUE_PAIRS[i][3],HUE_PAIRS[i][4],HUE_PAIRS[i][5],HUE_PAIRS[i][6],HUE_PAIRS[i][7]);
		}
		//--------------------------------------------------------------------------------------------------
	fprintf(tab,"\n");
		for(i=1;i<=48;i++)
		{
			    //int A[49],Aa[77];
				MAX_TEL[i].T = _read_data(inp);
				MAX_TEL[i].V = s = _read_data(inp);
				MAX_TEL[i].Y = 0.01 * s * s;
				fprintf(tab,"\n\t\ MAX  A =%3ld\tT =%5.1lf\tV =%5.1lf\tY =%5.1lf",
					A[i], MAX_TEL[i].T, MAX_TEL[i].V, MAX_TEL[i].Y);
		}

	fprintf(tab,"\n");
		//TELITETTSEG_GRAY[8][6][6];

		for(i=1; i <=7; i++)
		for(j=1; j <=5; j++)
		for(k=1; k <=5; k++)
			TELITETTSEG_GRAY[i][j][k] = 0;

		for(i=1; i <=7; i++)
		{

			fprintf(tab,"\n");
			i1 = 10*i;
			for(j=4; j >=2; j--)
			{
				fprintf(tab,"\n");
				for(k = 1; k <=5;  k++)
				{
					_read_data(inp);
					_read_data(inp);
					
					s = _read_data(inp);
					//if(k==1) TELITETTSEG_LEPCSO[i][j] = s;
					//helyette a max solid color telitettseg 5-ode az egyseg

					TELITETTSEG_GRAY[i][k][j] = _read_data(inp);
					fprintf(tab,"\n\t\tA =%3ld\tV = %5.2lf\t rel szurke = %5.2lf",
					10*i,5. + j*20., TELITETTSEG_GRAY[i][k][j] /*ATV sorrend*/);


				}//k
			}//j
		}//i

				VIL_GRAY[5] = 100.;  //Y = 100
				VIL_GRAY[4] =  85.;  //Y = 72.25
				VIL_GRAY[3] =  65.;  //Y = 42.25
				VIL_GRAY[2] =  45.;  //Y = 20.25
				VIL_GRAY[1] =   0.;  //Y = 0


		fclose(inp);
		fclose(tab);

		return(0);
}



