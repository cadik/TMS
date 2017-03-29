//magassagi szinskala vagy hamis szinezesi skala

//Lab		Selected_Colors[256];
//long	Selected_Heights[256];
//double	Relative_Delta[255];
//Lin_rgb _rgb[2600];
//meterenkent a magassagi rgb linearis szinek

void READ_MAGASSAGI_SZINEK(long * Num_of_Selected_Colors,
						   long * Num_of_Height_levels,
						   long * Num_of_Color_Sequence)
{
//magassagi intervallumok szama

	FILE * adat;
	long i, k;
	double X, Y, Z, R, G, B;

	adat = fopen("magassagi_szin.txt","rt");
	//adat = fopen("csont_hus.txt","rt");

	k=0;
	while(1)
	{
		AA[k]	=	_read_data(adat);
		TT[k]	=	_read_data(adat);
		VV[k]	=	_read_data(adat);

		if(AA[k] < 0) break;

		ATV_XYZ(AA[k], TT[k], VV[k], &X, &Y, &Z);
		XYZ_RGB709(X, Y, Z, &R, &G, &B);
		if(R<0 || R>1 ||
		   G<0 || G>1 ||
		   B<0 || B>1) 
		{
		   printf("\n\t%3ld. color is not displayble, RGB = (%6.3lf,%6.3lf,%6.3lf)",
					k,R,G,B);
		   getch();
		}

		XYZ_Hunter_Lab(	X,	Y,	Z, 
						&Selected_Colors[k].L,
						&Selected_Colors[k].a,
						&Selected_Colors[k].b);

		printf("\n A = %5.1f  T = %5.1f  V = %5.1f",
			       AA[k],TT[k],VV[k]);
		printf("\n\t %2ld  Lab   %lf  %lf  %lf",k,
						Selected_Colors[k].L,
						Selected_Colors[k].a,
						Selected_Colors[k].b);
		k++;
		*Num_of_Selected_Colors = k;
		//1-gyel kisebb a szinek szamanal
	}

	fclose(adat);

	//adat = fopen("_Selected_Heights.txt","rt");
	adat = fopen("csont_hus_density.txt","rt");
	*Num_of_Height_levels = _read_data(adat);
	//magassagi szintek 
	printf("\n  levels = %ld",*Num_of_Height_levels);

	for(k=0;k<(*Num_of_Height_levels);k++)
	{
		Selected_Heights[k] = (long)_read_data(adat);
		if(k==0 || k==(* Num_of_Height_levels)-1)
			printf("\n  %ld level  = %ld",k,Selected_Heights[k]);
	}

	fclose(adat);

	for(k=1;k<(*Num_of_Height_levels);k++)
	{
		delta[k]=1./(Selected_Heights[k] - Selected_Heights[k-1]);
		
		for(i=Selected_Heights[k-1];i<Selected_Heights[k];i++)
		{
			Height_delta[i-Selected_Heights[0]] = delta[k];
		//	printf("\n index = %ld delta = %lf",
		//		i-Selected_Heights[0],Height_delta[i-Selected_Heights[0]]);
		}
	}
	//a magassagi lepcsok relativ sziningerkulonbsege
	//Height_delta[i] a rel col diff [i,i+1] magassagintervallumban
	//i=> Selected_Heights[0],
	//    Selected_Heights[*Num_of_Height_levels-1]-1)

	*Num_of_Color_Sequence = 
		Selected_Heights[*Num_of_Height_levels-1]-
					   Selected_Heights[0];
	printf("\n  Num_of_Color_Sequence = %ld",
		*Num_of_Color_Sequence);
	
	printf("\n\teddig OK");
	getch();

}

void pseudo_colors(void)
{
	FILE * out;
	double dE, R, G, B, X, Y, Z, L, a, b, Xn, Yn, Zn;
	long i1, i2, i, j, k, intervallumok_szama;
	long Num_of_Selected_Colors, Num_of_Height_levels,
		 Num_of_Color_Sequence;

	Xn = D65_X;
	Yn = 100.;
	Zn = D65_Z; 


	READ_MAGASSAGI_SZINEK(&Num_of_Selected_Colors,
						  &Num_of_Height_levels,
						  &Num_of_Color_Sequence);

	for(k=1; k <= 7; k++)
	{
		dL[k] = (Selected_Colors[k].L - Selected_Colors[k-1].L) / delta[k];
		da[k] = (Selected_Colors[k].a - Selected_Colors[k-1].a) / delta[k];
		db[k] = (Selected_Colors[k].b - Selected_Colors[k-1].b) / delta[k];

		//dE = sqrt(dL[k]*dL[k]+da[k]*da[k]+db[k]*db[k]);
		printf("\n delta L = %lf  a = %lf   b = %lf", //  E = %lf",
			dL[k],da[k],db[k]); //,dE);
	}

	for(k=1; k <= Num_of_Color_Sequence; k++)
	{
	
		i1 = (long)Selected_Heights[k-1];
		i2 = (long)Selected_Heights[k];

		for(i=i1;i<i2;i++)
		{
			j = i - i1;
			L = Selected_Colors[k-1].L + j*dL[k];
			a = Selected_Colors[k-1].a + j*da[k];
			b = Selected_Colors[k-1].b + j*db[k];

			Lab_XYZ(L, a, b, Xn, Yn, Zn, &X, &Y, &Z);
			XYZ_RGB709(X, Y, Z, &R, &G, &B);
			RGBlinear_RGBdisplay(R, G, B, &Rps[i], &Gps[i], &Bps[i]);
		}

	}

	out = fopen("pseudo_color_sRGB.dat","wb");

	for(i=0;i<256;i++)
	{
		fwrite(&Rps[j],1,1,out);
		fwrite(&Gps[j],1,1,out);
		fwrite(&Bps[j],1,1,out);

		if(i%10==0)
		{
			printf("\n\tRGB  i = %4ld j = %4ld  %lf  %3lf  %3lf",
					i,j,Rps[j],Gps[j],Bps[j]);
			getch();
		}

	}

	fclose(out);

}

void pseudo_colors_2(void)
{
	FILE * out;
	double X,Y,Z,R,G,B;
	unsigned char rr, gg, bb;
	long i;
	double At,Bt, Av,Bv;

	out = fopen("Novak23.dat","wb");

	ATV_XYZ(65., 13., 42., &X, &Y, &Z);
	XYZ_RGB709(X, Y, Z, &R, &G, &B);
	RGBlinear_RGBdisplay(R, G, B, &rr, &gg, &bb);//0

	fwrite(&rr,1,1,out);
	fwrite(&gg,1,1,out);
	fwrite(&bb,1,1,out);

	for(i=1;i<=254;i++)
	{
	At = (8.-32.)/253.;
	Bt = 32. - At;

	Av = (74.-30.)/253.;
	Bv = 30. - Av;

	ATV_XYZ(23., At*(double)i+Bt, Av*(double)i+Bv, &X, &Y, &Z);
	XYZ_RGB709(X, Y, Z, &R, &G, &B);
	RGBlinear_RGBdisplay(R, G, B, &rr, &gg, &bb);//i

	fwrite(&rr,1,1,out);
	fwrite(&gg,1,1,out);
	fwrite(&bb,1,1,out);

	}

	ATV_XYZ(65., 22., 92., &X, &Y, &Z);
	XYZ_RGB709(X, Y, Z, &R, &G, &B);
	RGBlinear_RGBdisplay(R, G, B, &rr, &gg, &bb);//255

	fwrite(&rr,1,1,out);
	fwrite(&gg,1,1,out);
	fwrite(&bb,1,1,out);

	fclose(out);

}

long elso_index(long next_index,
				long Num_of_Selected_Colors,
				double del)
{
//az adott pontbol az elso szakaszveg az aktual_index lenne
//addig megyunk, amig del nagysagu ugras nem lehetseges
//return = -1, akkor elfogy a torottvonal
//return = a legkozelebbi jo szin-szakasz vegpontjanak indexe

	long k;
	double s2, dL, da, db;

	k = next_index;
	while(1)
	{
		dL = Selected_Colors[k].L - actual_point.L;
		da = Selected_Colors[k].a - actual_point.a;
		db = Selected_Colors[k].b - actual_point.b;

		s2 = dL*dL + da*da + db*db;
		if(s2 >= del*del) return(k);
		k++;
		if(k > Num_of_Selected_Colors) return(-1);
	}

}

void new_height_interval(long color_interval_index,
						 long height_interval_index,
						 long vegpont_index,
						 double del)
{
	double s, k2, v2, kv;

	if(color_interval_index == vegpont_index)
	{
	//azonos color szakaszon araszolunk, 
	//a kezdo es vegpont azonos color szakaszon van

	actual_point.L += del * unit_vector.L;
	Color_Lab_sequence[height_interval_index].L = actual_point.L;

	actual_point.a += del * unit_vector.a;		
	Color_Lab_sequence[height_interval_index].a = actual_point.a;

	actual_point.b += del * unit_vector.b;
	Color_Lab_sequence[height_interval_index].b = actual_point.b;
	}
	else
	{
	//itt tores utani pontot keresunk, 
	//egyes eseteben tobb color szakasszal kesobbi helyen

		k_vector.L = Selected_Colors[color_interval_index].L 
				   - actual_point.L;
		k_vector.a = Selected_Colors[color_interval_index].a 
				   - actual_point.a;
		k_vector.b = Selected_Colors[color_interval_index].b
				   - actual_point.b;	

		unit_vector.L = Selected_Colors[vegpont_index].L 
					  - Selected_Colors[vegpont_index-1].L;
		unit_vector.a = Selected_Colors[vegpont_index].a 
					  -	Selected_Colors[vegpont_index-1].a;
		unit_vector.b = Selected_Colors[vegpont_index].b
					  -	Selected_Colors[vegpont_index-1].b;	

		v2  =   unit_vector.L*unit_vector.L + 
				unit_vector.a*unit_vector.a + 
				unit_vector.b*unit_vector.b;

		s = 1./sqrt(v2);

		unit_vector.L *= s;
		unit_vector.a *= s;
		unit_vector.b *= s;

		k2  =	k_vector.L*k_vector.L + 
				k_vector.a*k_vector.a + 
				k_vector.b*k_vector.b;

		kv  =   k_vector.L*unit_vector.L + 
				k_vector.a*unit_vector.a + 
				k_vector.b*unit_vector.b;

		s = - kv + sqrt(kv * kv + del*del - k2);

		actual_point.L = Selected_Colors[vegpont_index-1].L
					   + s * unit_vector.L;
		Color_Lab_sequence[height_interval_index].L = actual_point.L;

		actual_point.a = Selected_Colors[vegpont_index-1].a
			           + s * unit_vector.a;
	    Color_Lab_sequence[height_interval_index].a = actual_point.a;

		actual_point.b = Selected_Colors[vegpont_index-1].b
			           + s * unit_vector.b;
		Color_Lab_sequence[height_interval_index].b = actual_point.b;
		//ezt be kell irni a vegso szinsorba
	}
	
//	printf("\n vegpont = %3ld  Lab = %7.3lf  %7.3lf  %7.3lf",
//		 vegpont_index,actual_point.L,actual_point.a,actual_point.b);

//szin tombbe beirjuk az aktual-t

}

double remain_color_difference(long vegpont_index,
							   long Num_of_Selected_Colors)
{
//ekkor a legnagyobb magassagig eljutottunk,
//de nem jartuk be a szin torottvonalat
//return a maradek torottvonal szakasz ossz-hossza
//ami col diff ertekek osszege
	long k;
	double s, sL, sa, sb, t;

	sL = Selected_Colors[vegpont_index].L - actual_point.L;
	sa = Selected_Colors[vegpont_index].a - actual_point.a;
	sb = Selected_Colors[vegpont_index].b - actual_point.b;

	s = sqrt(sL*sL + sa*sa +sb*sb);
	if(vegpont_index == Num_of_Selected_Colors) return(s);

	for(k=vegpont_index+1;k<=Num_of_Selected_Colors;k++)
	{
		unit_vector.L = Selected_Colors[k].L 
					  - Selected_Colors[k-1].L;
		unit_vector.a = Selected_Colors[k].a
			          - Selected_Colors[k-1].a;
		unit_vector.b = Selected_Colors[k].b
			          - Selected_Colors[k-1].b;

		t=sqrt(unit_vector.L * unit_vector.L + 
			   unit_vector.a * unit_vector.a + 
			   unit_vector.b * unit_vector.b);
		s += t;
	}

	return(s);

}


double tortvonal_bejarasa_adott_konstanssal(
				double konstans,
				long Num_of_Selected_Colors,
				long Num_of_Height_levels,
				long Num_of_Color_Sequence)
{
//return = -1, akkor tul nagy a konstans
//return = utolso magassagpont es az utolso szinpont tavolsaga
//ha return < EPS, akkor elfogadhato a bejaras

	long vegpont_index, //max_height_value,
		 color_interval_index, height_index;
	double s, del;
	//unit_vector, actual_point

	color_interval_index  = 1;
	height_index = 0;
    //max_height_value = Selected_Heights[Num_of_Height_levels-1]-
	//				   Selected_Heights[0];

	printf("\n  REL  min_H =%ld  MAX=%ld",
				0, Num_of_Color_Sequence);

	Color_Lab_sequence[0].L = 
		    actual_point.L = Selected_Colors[0].L;
	Color_Lab_sequence[0].a = 
			actual_point.a = Selected_Colors[0].a;
	Color_Lab_sequence[0].b = 
			actual_point.b = Selected_Colors[0].b;

	while(1)
	{
	//szakaszonkent jarjuk be a szinpontok sorozatat
	
		unit_vector.L = Selected_Colors[color_interval_index].L 
					  - Selected_Colors[color_interval_index-1].L;
		unit_vector.a = Selected_Colors[color_interval_index].a
			          - Selected_Colors[color_interval_index-1].a;
		unit_vector.b = Selected_Colors[color_interval_index].b
			          - Selected_Colors[color_interval_index-1].b;

		s = 1. / sqrt(unit_vector.L * unit_vector.L + 
					  unit_vector.a * unit_vector.a + 
					  unit_vector.b * unit_vector.b);
		
		unit_vector.L *= s;
		unit_vector.a *= s;
		unit_vector.b *= s;
		//a bejaras iranyaba mutato normal egysegvektor

		del = konstans * Height_delta[height_index];

		if((vegpont_index = 
		    elso_index( color_interval_index, 
						Num_of_Selected_Colors, 
						del ) ) < 0)
			return(-1.);

		//ekkor tul nagy volt a felvett konstans
		//vegigaraszoltuk az egesz torottvonalat es elfogyott
		//ha itt vagyunk, akkor delta hosszu szakasz
		//vegpontjat keressuk a [vegpont_index-1, vegpont_index]
		//Selected_Colors szakaszon az actual_pointbol indulva

		new_height_interval(color_interval_index, 
							height_index,
							vegpont_index,
							del);

		printf("\n\tH_index = %ld  vegpont_index = %ld",
			height_index,vegpont_index);

		color_interval_index = vegpont_index;
		if(height_index < Num_of_Color_Sequence)
				height_index++;
		else
			return(remain_color_difference(
			vegpont_index, Num_of_Selected_Colors));
	}
}
/*
void magassagi_rgb_sorozat_(void)
{
	FILE * out;
	double dE, R, G, B, X, Y, Z, L, a, b, Xn, Yn, Zn;
	long i1, i2, i, j, k, intervallumok_szama;

	Xn = D65_X;
	Yn = 100.;
	Zn = D65_Z;

	READ_MAGASSAGI_SZINEK(&intervallumok_szama);

	for(k=1; k <= intervallumok_szama; k++)
	{
		dm[k] = magassag[k] - magassag[k-1];
		dL[k] = (LL[k] - LL[k-1]) / dm[k];
		da[k] = (aa[k] - aa[k-1]) / dm[k];
		db[k] = (bb[k] - bb[k-1]) / dm[k];

		dE = sqrt(dL[k]*dL[k]+da[k]*da[k]+db[k]*db[k]);
		printf("\n delta L = %lf  a = %lf   b = %lf  E = %lf",
			dL[k],da[k],db[k],dE);
	}

	for(k=1; k <= intervallumok_szama; k++)
	{
	
		i1 = (long)magassag[k-1];
		i2 = (long)magassag[k];
		for(i=i1;i<i2;i++)
		{
			j = i - i1;
			L = LL[k-1] + j*dL[k];
			a = aa[k-1] + j*da[k];
			b = bb[k-1] + j*db[k];

			Lab_XYZ(L, a, b, Xn, Yn, Zn, &X, &Y, &Z);
			XYZ_RGB709(X, Y, Z, &R, &G, &B);
			r_[i] = R;
			g_[i] = G;
			b_[i] = B;
		}

	}

	ATV_XYZ(52., 7., 35., &X, &Y, &Z); //háttér színe
	XYZ_RGB709(X, Y, Z, &R, &G, &B);
			r_[0] = R;
			g_[0] = G;
			b_[0] = B;
		
	ATV_XYZ(56., 30., 70., &X, &Y, &Z); //víz színe
	XYZ_RGB709(X, Y, Z, &R, &G, &B);
			r_[1] = R;
			g_[1] = G;
			b_[1] = B;

	ATV_XYZ(51.5, 11., 100., &X, &Y, &Z); //egbolt színe
	XYZ_RGB709(X, Y, Z, &R, &G, &B);
			r_[2] = R;
			g_[2] = G;
			b_[2] = B;
	printf("\n RGB = %e  %e  %e",R,G,B);

	ATV_XYZ(20.124, 11., 100., &X, &Y, &Z); //napfeny színe
	XYZ_RGB709(X, Y, Z, &R, &G, &B);
			r_[3] = R;
			g_[3] = G;
			b_[3] = B;
	printf("\n RGB = %e  %e  %e",R,G,B);
	
	ATV_XYZ(22., 15., 60., &X, &Y, &Z); //szintvonalak színe
	XYZ_RGB709(X, Y, Z, &R, &G, &B);
			r_[4] = R;
			g_[4] = G;
			b_[4] = B;
	printf("\n RGB = %e  %e  %e",R,G,B);
	

	out = fopen("magassagi_rgb.dat","wb");

	for(i=0;i<2000;i++)
	{
		fwrite(&r_[j],1,sizeof(double),out);
		fwrite(&g_[j],1,sizeof(double),out);
		fwrite(&b_[j],1,sizeof(double),out);

		if(i%100==0)
		{
			printf("\n\tRGB  i = %4ld j = %4ld  %lf  %3lf  %3lf",
					i,j,r_[j],g_[j],b_[j]);
			getch();
		}
	}

	fclose(out);

}
*/
void Bartos_(void)
{

	FILE * adat , * out;
	long i;
	double a, t, v, X, Y, Z, R, G, B;
	unsigned char rr, gg, bb;

	adat = fopen("Bartos.txt","rt");
	out  = fopen("Bartos_sRGB.txt","wt");

	for(i=0;i<230;i++)
	{
		a	=	_read_data(adat);
		t	=	_read_data(adat);
		v	=	_read_data(adat);

		ATV_XYZ(a, t, v, &X, &Y, &Z);
		XYZ_RGB709(X, Y, Z, &R, &G, &B);
		RGBlinear_RGBdisplay(R, G, B, &rr, &gg, &bb);

		fprintf(out,"\n\t%u\t%u\t%u\t%3ld\t%3ld\t%3ld",
			rr,gg,bb,(long)a,(long)t,(long)v);

	}

	fclose(out);
	fclose(adat);


}