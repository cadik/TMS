// TMOCmd.cpp: implementation of the TMOCmd class.
//
//////////////////////////////////////////////////////////////////////
#include "TMOCmd.h"

#ifdef LINUX
 #include "../TMOLinux/TMOLinux.h"
#endif
#ifndef LINUX
 #include "../TMOW32/TMOW32.h"
#endif


#include <stdlib.h>
#include <wchar.h>
#include <iostream>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOCmd::TMOCmd()
{

}

TMOCmd::~TMOCmd()
{

}

int TMOCmd::main(int argc, char *argv[])
{
	int num_libraries=0, i, j, pcount, temp, opindex;
	std::string tempStr;
	int parameter_id, picture_count = 0;
	double d;
	TMOImage input;
	TMOImage output;
	TMO** op;
	TMOParameter **params;
	wchar_t buffer[256];
	char buffer1[20];
	wchar_t* pLib[256];

	op = new TMO*[256];
	num_libraries = EnumLibraries(pLib, 256);
	iOperatorCount = 0;
	for (j = 0; j < num_libraries; j++)
	{
		iOperatorCount += OpenLibrary(pLib[j], &op[iOperatorCount]);
	}

	wprintf (L"TMO - Tone Mapping Operator application\n\n");

	if (argc < 2) 
	{
		Help(op);
		return 3;
	}
	for (opindex = 0; opindex < iOperatorCount; opindex++)
	{
		mbstowcs(buffer,argv[1],255);
		if (!wcscmp(buffer, op[opindex]->GetName()))
		{
			wprintf (L"%ls\n\n", op[opindex]->GetDescription());
			pcount = op[opindex]->GetParameterCount();
			params = new TMOParameter*[pcount];
			op[opindex]->EnumParameters(params);
			break;
		}
	}
	if (opindex ==  iOperatorCount) 
	{
		fwprintf( stderr, L"Wrong tone mapping method : %s\n",buffer);
		Help(op);
		return 2;
	}

	for (i = 2; i < argc; i++)
	{
		if (argv[i][0] == '-')
		{
			parameter_id = 0;
			if (strlen(argv[i])==1) 
			{
				fprintf( stderr, "Wrong parameter : -\n");
				parameter_id = -1;
			}
			mbstowcs(buffer,&argv[i][1],255);
			for (j = 0; j < pcount; j++)
			{
				if (!wcscmp(buffer, params[j]->GetName()))
				{
					if (params[j]->Is(TMO_BOOL))
					{
						if (argv[i+1][0] == '-')
						{
							*params[j] = true;
						}
						else
						{
							i++;
							if ((strcmp(argv[i], "On")==0)||(strcmp(argv[i], "ON")==0)||
								(strcmp(argv[i], "on")==0)||(strcmp(argv[i], "1")==0))
								*params[j] = true;
							else if ((strcmp(argv[i], "Off")==0)||(strcmp(argv[i], "OFF")==0)||
								(strcmp(argv[i], "off")==0)||(strcmp(argv[i], "0")==0))
								*params[j] = false;
							else 
							{
								*params[j] = true;
								i--;
							}
						}
						if (params[j]->GetBool()) wprintf (L"- %s : ON\n",params[j]->GetDescription());
						else wprintf (L"- %s : OFF\n",params[j]->GetDescription());
					}
					if (params[j]->Is(TMO_DOUBLE))
					{
						i++;
						d = strtod(argv[i],0);
						wprintf (L"- %ls : %f\n",params[j]->GetDescription(),d);
						*params[j] = d;
					}
					if (params[j]->Is(TMO_INT))
					{
						i++;
						temp = strtol(argv[i],0,10);
						wprintf (L"- %ls : %i\n",params[j]->GetDescription(),temp);
						*params[j] = temp;
					}
					if (params[j]->Is(TMO_STRING))
					{
						i++;						
						tempStr = argv[i];
						wprintf (L"- %ls : %s\n",params[j]->GetDescription(), tempStr.c_str());
						*params[j] = tempStr;
					}
				}
			}
		}
		else 
		{
			picture_count = 1;
			break;
		}
	}
	if (!picture_count) 
	{
		Help(op, opindex);
		return 1;
	}

	try
	{
		input.Open(argv[i]);
		input.Convert(TMO_RGB);
		//input.SaveWithSuffix("_HDRinput", TMO_RAW);
		//input.SaveWithSuffix("_HDRinput", TMO_EXR);

		output.New(input);
	
		op[opindex]->SetImages(input, output);
		op[opindex]->Transform();
		
		output.Convert(TMO_RGB);
		wcstombs(buffer1, op[opindex]->GetName(), 20);
		output.SaveWithSuffix(buffer1);
		//output.SaveWithSuffix(buffer1, TMO_RAW);
		output.SaveWithSuffix(buffer1, TMO_EXR_16);
	}
	catch(int a)
	{
		if (a == -1) fprintf(stderr,"File not found : %s\n", argv[i]);
		else fprintf (stderr,"Exception number %i.\n",a);
	}
	
	iOperatorCount = 0;
	for (j = 0; j < num_libraries; j++)
	{
		iOperatorCount += CloseLibrary(pLib[j], &op[iOperatorCount]);
	}

	delete[] op;
	return 0;
}

void TMOCmd::Help(TMO** op, int opindex)
{
	int i, j, count;
	TMOParameter** params;
	wprintf (L"\nUsage : \n\n");
	wprintf (L"    ./tmocmd method [-arg1 -arg2 ... -argn] filename\n\n\n");
	wprintf (L"Available methods :\n\n");
	for (j = 0; j < iOperatorCount; j++)
	{
		if ((opindex > 0)&&(opindex!=j)) continue;
		wprintf (L"    %ls ... %ls\n\n", op[j]->GetName(), op[j]->GetDescription());
		wprintf (L"        Parameters : \n\n");
		count = op[j]->GetParameterCount();
		params = new TMOParameter*[count];
		op[j]->EnumParameters(params);
		for (i = 0; i < count; i++)
		{
			if (params[i]->Is(TMO_BOOL))
				wprintf(L"            -%ls $b$ ... %ls\n", params[i]->GetName(), params[i]->GetDescription());
			if (params[i]->Is(TMO_DOUBLE))
				wprintf(L"            -%ls $d$ ... %ls\n", params[i]->GetName(), params[i]->GetDescription());
			if (params[i]->Is(TMO_INT))
				wprintf(L"            -%ls $i$ ... %ls\n", params[i]->GetName(), params[i]->GetDescription());
			if (params[i]->Is(TMO_STRING))
				wprintf(L"            -%ls $s$ ... %ls\n", params[i]->GetName(), params[i]->GetDescription());
		}
		wprintf (L"\n\n");
		delete[] params;
	}
	wprintf (L"\n    $b$ ... boolean value (On, Off, 0, 1)\n");
	wprintf (L"    $d$ ... floating-point value\n");
	wprintf (L"    $i$ ... integer value\n");
	wprintf (L"    $s$ ... string\n\n");
}

int TMOCmd::Bar(int part, int all)
{
	if (part==all) printf ("100%%\n");
	else printf ("%2i%%\b\b\b", 100*part/all);
	return 0;
}
