
//GLOBAL long _2D_uniform_matrix[256][256];

void dithering_preprocessing(void)
{

FILE * dit;

long x, y, i, j;

x = y = 256;

dit = fopen("dither_256_256.dat","rb");

for(i=0;i<x;i++)
	for(j=0;j<y;j++)
		fread(&_2D_uniform_matrix[i][j],4,1,dit);

fclose(dit);

}

long hanyadik_szin_van_a_kerdeses_negyzeten(
				long x_offset,
				long y_offset, //idofuggo random [0,200] offset ertekek
				long x_negyzet,
				long y_negyzet, //a kerdeses szinu negyzet koordinatai
				long szinfajtak_szama_a_szincsoportban) 
					//N, ekkor 0,1,...,(N-1) a indexu 
					//szinekbol valasztahatunk, ezek a szincsoport elemei
{
//return = a lehetseges 0,1,...,N-1 szinsorszamok egyike
	long d, csoport_meret, szin_sorszam;

	d = _2D_uniform_matrix[x_offset+x_negyzet][y_offset+y_negyzet];

	csoport_meret = 65536/szinfajtak_szama_a_szincsoportban;

	szin_sorszam = (long) (d / csoport_meret);

	if(szin_sorszam >= szinfajtak_szama_a_szincsoportban) 
	   szin_sorszam  = szinfajtak_szama_a_szincsoportban - 1;

	return(szin_sorszam);

}

//pelda: a (9,5) indexu negyzet szinet keressuk a 4 szinbol
//hanyadik_szin_van_a_kerdeses_negyzeten(0,0,9,5,4)
//itt offset az egyszeruseg kedveert nincs
//return-kent a 0,1,2 vagy 3 szin sorszamok egyiket kapjuk
