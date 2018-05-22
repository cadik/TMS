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

#ifndef LINUX
 #define OPENEXR_DLL
#else
 #define _stricmp strcasecmp
#endif


#include <stdlib.h>
#include <limits.h>
#include <wchar.h>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOCmd::TMOCmd()
{

}

TMOCmd::~TMOCmd()
{

}

void closeLibraries(int numLibraries, wchar_t** libs, TMO** operators){
	int iOperatorCount = 0;
	for (int j = 0; j < numLibraries; j++)
	{
		iOperatorCount += CloseLibrary(libs[j], &operators[iOperatorCount]);
	}
}

int argtoi(char *arg, int *value){
	char *ptr=NULL;
	long int l = strtol(arg,&ptr,10);
	if(*ptr!='\0' || l>INT_MAX || l<INT_MIN)
		return 1;
	
	*value = l;
	return 0;
}

int argtod(char *arg, double *value){
	char *ptr=NULL;
	*value = strtod(arg,&ptr);
	if(*ptr!='\0')
		return 1;
	return 0;
}

int TMOCmd::main(int argc, char *argv[])
{
	int num_libraries=0, i, j, pcount, temp, opindex;
	std::string tempStr;
	int picture_count = 0, pfound;
	double d;
	TMOImage input;
	TMOImage output;
	TMOVideo inVideo;
	TMOVideo outVideo;
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
		closeLibraries(num_libraries,pLib,op);
		delete[] op;
		return 3;
	}
	for (opindex = 0; opindex < iOperatorCount; opindex++)
	{
		mbstowcs(buffer,argv[1],255);
		if (!wcscmp(buffer, op[opindex]->GetName()))
		{
			pcount = op[opindex]->GetParameterCount();
			params = new TMOParameter*[pcount];
			op[opindex]->EnumParameters(params);
			break;
		}
	}
	if (opindex ==  iOperatorCount) 
	{
		Help(op);
		fwprintf(stderr, L"Wrong tone mapping method : %ls\n",buffer);
		closeLibraries(num_libraries,pLib,op);
		delete[] op;
		return 2;
	}

	wprintf (L"%ls\n\n", op[opindex]->GetDescription());
	if (argc == 2){
		Help(op, opindex);
		closeLibraries(num_libraries,pLib,op);
		delete[] op;
		return 0;
	} 
	
	for (i = 2; i < argc; i++)
	{
		if (argv[i][0] == '-')
		{
			pfound = 0;
			mbstowcs(buffer,&argv[i][1],255);
			for (j = 0; j < pcount; j++)
			{
				if (!wcscmp(buffer, params[j]->GetName()))
				{
					pfound = 1;
					if (params[j]->Is(TMO_BOOL))
					{
						if (i+1 == argc) break;
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
						if (params[j]->GetBool()) wprintf (L"- %ls : ON\n",params[j]->GetDescription());
						else wprintf (L"- %ls : OFF\n",params[j]->GetDescription());
					}
					
					if(i+1 == argc){
						fwprintf( stderr, L"Missing value for parameter : %s\nTry: '%s %ls' for more information\n",argv[i],argv[0],op[opindex]->GetName());
						closeLibraries(num_libraries,pLib,op);
						delete[] op;
						return 5;
					}
					
					if (params[j]->Is(TMO_DOUBLE))
					{
						i++;
						if(argtod(argv[i],&d)){
							fwprintf( stderr, L"Wrong floating-point format : %s, for parameter : %s\nTry: '%s %ls' for more information\n",argv[i],argv[i-1],argv[0],op[opindex]->GetName());
							closeLibraries(num_libraries,pLib,op);
							delete[] op;
							return 6;
						}
						wprintf (L"- %ls : %f\n",params[j]->GetDescription(),d);
						*params[j] = d;
					}
					if (params[j]->Is(TMO_INT))
					{
						i++;
						if(argtoi(argv[i],&temp)){
							fwprintf( stderr, L"Wrong integer format : %s, for parameter : %s\nTry: '%s %ls' for more information\n",argv[i],argv[i-1],argv[0],op[opindex]->GetName());
							closeLibraries(num_libraries,pLib,op);
							delete[] op;
							return 6;
						}
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
			if(!pfound){
				fwprintf( stderr, L"Wrong parameter : %s\nTry: '%s %ls' for more information\n",argv[i],argv[0],op[opindex]->GetName());
				closeLibraries(num_libraries,pLib,op);
				delete[] op;
				return 4;
			}
		}
		else 
		{
			picture_count = 1;
			if(i+1!=argc)
				fprintf( stderr, "\nWarning : Additional parameters after '%s' ignored\n\n",argv[i]);
			break;
		}
	}
	if (!picture_count) 
	{
		fwprintf( stderr, L"Missing file operand\nTry: '%s %ls' for more information\n",argv[0],op[opindex]->GetName());
		closeLibraries(num_libraries,pLib,op);
		delete[] op;
		return 1;
	}
	

	try
	{
	  struct stat path_stat;
	  stat(argv[i], &path_stat);
	  
	  
	  if (S_ISREG(path_stat.st_mode)==1) //it is a file
	  {
	    int length = strlen(argv[i]), ii;
	    char* tmp;
	    //if (tmp) delete[] tmp;
	    tmp = new char[length + 6];
	    strcpy(tmp, argv[i]);

	    if (!length) throw TMO_EFILE;
	    ii = length - 1;
	    while (ii)
	    {
		    if (tmp[ii] == '.') break;
		    ii--;
	    }
	    ii++;
	    if(_stricmp(&tmp[ii],"mp4") == 0 || _stricmp(&tmp[ii],"avi") == 0 || _stricmp(&tmp[ii],"wmv") == 0) ///if video that output will be video
	    {
	      inVideo.OpenVideo(argv[i]);
  
	      wcstombs(buffer1, op[opindex]->GetName(), 20);
	      inVideo.setNameOut(buffer1);
	      outVideo.createOutputVideo(inVideo);
	     
	 
	      if(TMOv* opVid = dynamic_cast<TMOv*>(op[opindex])) //if TransforVideo is implemented cast and call
	      {
		
		  opVid->SetVideos(inVideo,outVideo);
		  opVid->TransformVideo(); 
	      }
	      else   ///if trasform video is not implemented cast will fail and frame by frame will be called
	      {
		for(int j=0;j<inVideo.GetTotalNumberOfFrames();j++)
		{
		  inVideo.getTMOImageVideoFrame(inVideo.getVideoCaptureObject(),j,input);
		  
		  input.Convert(TMO_RGB);	
		  output.New(input);
	  
		  op[opindex]->SetImages(input, output);  
		  op[opindex]->Transform();
		  
		  outVideo.setTMOImageFrame(outVideo.getVideoWriterObject(),output);
		  output.Close();
		}
	      }
   
	    }
	    else
	    {    //output will be image
		  
		  input.Open(argv[i]);
		  input.Convert(TMO_RGB);

		  output.New(input);
	  
		  op[opindex]->SetImages(input, output);  
		  op[opindex]->Transform();
		  
		  output.Convert(TMO_RGB);
		  wcstombs(buffer1, op[opindex]->GetName(), 20);
		  output.SaveWithSuffix(buffer1);
		  output.SaveWithSuffix(buffer1, TMO_EXR_16);
	     }
	  
	  }
	  else if (S_ISDIR(path_stat.st_mode)==1) // is a directory, the directory will be searched anad a video will be output
	  {
	    DIR* pDir = opendir ( argv[i]);
	    if ( !pDir) throw TMO_EFILE;
	    dirent* pEntry;
	    std::vector<std::string> fileNames;
	    std::string file;
	    std::string fileType;
	    
	    while ( pEntry = readdir ( pDir)) 
	    {
		file = pEntry->d_name;  //get name of file in dir
		    
		if(file.size()>4) //name of file must be bigger because for exmaple a.jpg
		{
		    fileType = file.substr( file.length() - 3 );
		    if (fileType == "hdr" || fileType == "pic" || fileType == "pfm" || fileType == "exr" ||
			fileType == "hdrraw" || fileType == "raw" || fileType == "tif" || fileType == "jpg" ||
			fileType == "jpeg" || fileType == "png" || fileType == "ppm")
		    {
			file =  "/" + file;
			file = argv[i] + file;      ///if image than add path
			fileNames.push_back(file);
			
		    }
		}    
	    }	     
		
	    for (int  it = 0 ; it <fileNames.size(); ++it)  ///every image in dir is converted
	    {
	      input.Open(fileNames[it].c_str());  
	      input.Convert(TMO_RGB);
	      output.New(input);
 
	      op[opindex]->SetImages(input, output);
	      op[opindex]->Transform();    ///start conversion
	      
	      output.Convert(TMO_RGB);
	      wcstombs(buffer1, op[opindex]->GetName(), 20);
	      output.SaveWithSuffix(buffer1);
		
	      output.SaveWithSuffix(buffer1, TMO_EXR_16);
	    
	    }
	  }
	  else throw TMO_EFILE; ///if it isnt anything nor a file nor directory
	  
	}
	catch(TMOIMAGE_API a)
	{
		if (a == TMO_EFILE) fprintf(stderr,"File not found : %s\n", argv[i]);
		else fprintf (stderr,"Exception number %i.\n",a);
		closeLibraries(num_libraries,pLib,op);
		delete[] op;
		return 7;
	}
	
	closeLibraries(num_libraries,pLib,op);
	delete[] op;
	return 0;
}

void TMOCmd::Help(TMO** op, int opindex)
{
	int i, j, count;
	TMOParameter** params;
	wprintf (L"Usage : \n");
	wprintf (L"    ./tmocmd method [-arg1 -arg2 ... -argn] filename\n\n");
	wprintf (L"Available methods :\n");
	for (j = 0; j < iOperatorCount; j++)
	{
		if ((opindex >= 0)&&(opindex!=j)) continue;
		wprintf (L"    %ls ... %ls\n", op[j]->GetName(), op[j]->GetDescription());
		wprintf (L"        Parameters : \n");
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
		wprintf (L"\n");
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
