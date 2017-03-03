//kiserlet_17.h

void kiserlet_17_plus_2_2(void)
{

long sorhossz[] = 
{0, 8, 10, 8, 9, 9, 10, 10, 8, 8, 10, 10, 10, 10, 9, 8, 10, 8};

long i,k,n;
short j;
char c;
double	 X,Y,Z,s,
		_A,T,V,wh,bl,col,lam,fi;

FILE *inp, *tabl;

tabl= fopen("UNIFORM_SERIES","wt");

inp = fopen("Kor_Y_60.txt","rt");

		V = 10. * sqrt(60.);

		for(i=1;i<=48;i++)
		{
		_read_data(inp);
		_A = A[i];
		T = _read_data(inp);
		ATV_XYZ(_A, T, V, &X, &Y, &Z);
		printf("\n\tA = %8.2lf  fi = %8.2lf   T=%8.2lf  V = %8.2lf",
			_A,Coloroid_fi_D65[i],T,V);
		fprintf(tabl,"\n\t%9.4lf\t%9.4lf\t%9.4lf\t%9.4lf\t%6.3lf\t%6.3lf\t%6.3lf",
			Coloroid_fi_D65[i],_A,T,V,X,Y,Z);
		}
		fprintf(tabl,"\n%9.4lf",-1.);
	printf("\n ---------------------------------------------------");
	
fclose(inp);

inp = fopen("Kor_Y_70.txt","rt");

		V = 10. * sqrt(70.);

		for(i=1;i<=48;i++)
		{
		_read_data(inp);
		_A = A[i];
		T = _read_data(inp);
		ATV_XYZ(_A, T, V, &X, &Y, &Z);
		printf("\n\tA = %8.2lf  fi = %8.2lf   T=%8.2lf  V = %8.2lf",
			_A,Coloroid_fi_D65[i],T,V);
		fprintf(tabl,"\n\t%9.4lf\t%9.4lf\t%9.4lf\t%9.4lf\t%6.3lf\t%6.3lf\t%6.3lf",
			Coloroid_fi_D65[i],_A,T,V,X,Y,Z);
		}
		fprintf(tabl,"\n%9.4lf",-1.);
	printf("\n ---------------------------------------------------");
	
fclose(inp);

inp = fopen("kettoszinsor_atvxyz.txt","rt");

		for(i=1;i<=30;i++)
		{
		_read_data(inp);
		_A = _read_data(inp);
		T = _read_data(inp);
		V = _read_data(inp);
		//ATV_XYZ(_A, T, V, &X, &Y, &Z);
		X = _read_data(inp);
		Y = _read_data(inp);
		Z = _read_data(inp);

		XYZ_ATV(X,Y,Z,&_A,&T,&V,&wh,&bl,&col,&lam,&fi);
		printf("\n\tA = %8.2lf  fi = %8.2lf   T=%8.2lf  V = %8.2lf",
				_A,fi,T,V);
	
		fprintf(tabl,"\n\t%9.4lf\t%9.4lf\t%9.4lf\t%9.4lf\t%6.3lf\t%6.3lf\t%6.3lf",
			Coloroid_fi_D65[i],_A,T,V,X,Y,Z);
		}
		fprintf(tabl,"\n%9.4lf",-1.);
	printf("\n ---------------------------------------------------");

		for(i=1;i<=40;i++)
		{
		_read_data(inp);
		_A = _read_data(inp);
		T = _read_data(inp);
		V = _read_data(inp);
		//ATV_XYZ(_A, T, V, &X, &Y, &Z);
		X = _read_data(inp);
		Y = _read_data(inp);
		Z = _read_data(inp);

		XYZ_ATV(X,Y,Z,&_A,&T,&V,&wh,&bl,&col,&lam,&fi);
		printf("\n\tA = %8.2lf  fi = %8.2lf   T=%8.2lf  V = %8.2lf",
				_A,fi,T,V);
	
		fprintf(tabl,"\n\t%9.4lf\t%9.4lf\t%9.4lf\t%9.4lf\t%6.3lf\t%6.3lf\t%6.3lf",
			Coloroid_fi_D65[i],_A,T,V,X,Y,Z);
		}
		fprintf(tabl,"\n%9.4lf",-1.);
	printf("\n ---------------------------------------------------");

fclose(inp);


inp = fopen("mintak.txt","rt");

for(i=1;i<=17;i++)
{
	for(n=1;n<=sorhossz[i];n++)
	Xave[n] = Yave[n] = Zave[n] = 0;

	for(k=1;k<=4;k++)
	{
		for(n=1;n<=sorhossz[i];n++)
		{
			_read_data(inp);
			X = _read_data(inp);
			Y = _read_data(inp);
			Z = _read_data(inp);
			if(k==2 || k==3) s = 1;
			if(k==4) s = 3;
				if(k>1)
				{
					Xave[n] += s * X;
					Yave[n] += s * Y;
					Zave[n] += s * Z;
				}
			//k=1 gyenge minõségû
			//k=2 és k=3 elég jó 
			//k = 4 : a legjobb (Nemcsics fele)
			//	sorozatok, szinte etalon ertekuek
		}//end_n

	}//end_k
		
	for(n=1;n<=sorhossz[i];n++)
	{
		X = Xave[n] / 5;
		Y = Yave[n] / 5;
		Z = Zave[n] / 5;
		//súlyozott átlagok

		XYZ_ATV(X,Y,Z,
			&_A,&T,&V,&wh,&bl,&col,&lam,&fi);
		printf("\n\tA = %8.2lf  fi = %8.2lf   T=%8.2lf  V = %8.2lf",
				_A,fi,T,V);
		fprintf(tabl,"\n\t%9.4lf\t%9.4lf\t%9.4lf\t%9.4lf\t%6.3lf\t%6.3lf\t%6.3lf",
			fi,_A,T,V,X,Y,Z);
	}

	fprintf(tabl,"\n%9.4lf",-1.);
	printf("\n ---------------------------------------------------");
	
}//end_i_17


fclose(tabl);

}

void kiserlet_plus(void)
{

long i,k,n;
short j;
char c;
double	 X,Y,Z,s,
		_A,T,V,wh,bl,col,lam,fi;

FILE *inp, *tabl;

tabl= fopen("kettoszinsor","wt");

inp = fopen("kettoszinsor_atvxyz.txt","rt");

		for(i=1;i<=30;i++)
		{
		_read_data(inp);
		_A = _read_data(inp);
		T = _read_data(inp);
		V = _read_data(inp);
		//ATV_XYZ(_A, T, V, &X, &Y, &Z);
		X = _read_data(inp);
		Y = _read_data(inp);
		Z = _read_data(inp);

		XYZ_ATV(X,Y,Z,&_A,&T,&V,&wh,&bl,&col,&lam,&fi);
		printf("\n\tA = %8.2lf  fi = %8.2lf   T=%8.2lf  V = %8.2lf",
				_A,fi,T,V);
	
		fprintf(tabl,"\n\t%9.4lf\t%9.4lf\t%9.4lf\t%9.4lf\t%6.3lf\t%6.3lf\t%6.3lf",
			Coloroid_fi_D65[i],_A,T,V,X,Y,Z);
		}
		fprintf(tabl,"\n%9.4lf",-1.);
	printf("\n ---------------------------------------------------");

		for(i=1;i<=40;i++)
		{
		_read_data(inp);
		_A = _read_data(inp);
		T = _read_data(inp);
		V = _read_data(inp);
		//ATV_XYZ(_A, T, V, &X, &Y, &Z);
		X = _read_data(inp);
		Y = _read_data(inp);
		Z = _read_data(inp);

		XYZ_ATV(X,Y,Z,&_A,&T,&V,&wh,&bl,&col,&lam,&fi);
		printf("\n\tA = %8.2lf  fi = %8.2lf   T=%8.2lf  V = %8.2lf",
				_A,fi,T,V);
	
		fprintf(tabl,"\n\t%9.4lf\t%9.4lf\t%9.4lf\t%9.4lf\t%6.3lf\t%6.3lf\t%6.3lf",
			Coloroid_fi_D65[i],_A,T,V,X,Y,Z);
		}
		fprintf(tabl,"\n%9.4lf",-1.);

fclose(inp);

fclose(tabl);

}