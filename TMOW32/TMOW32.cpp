// TMOW32.cpp : Defines the entry point for the DLL application.
//
//#include "../TMOPlugin/TmoPlugin.h"
#include "../tmolib/TMO.h"
#include "TMOW32.h"
#include <map>
#include <string>

using namespace std;

map<wstring,HMODULE> mapModules;

#define MAX_OPERATORS 256

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    return TRUE;
}

int TMOW32_API EnumLibraries(wchar_t** filenames, int maxcount)
{
#ifdef WIN_NT
	WIN32_FIND_DATAW FileData; 
#else 
	WIN32_FIND_DATAA FileData; 
#endif
	HANDLE hSearch; 
	BOOL fFinished = FALSE; 
	int i = 0;


#ifdef WIN_NT
	hSearch = FindFirstFileW(L"*.tms", &FileData); 
#else 
	hSearch = FindFirstFileA("*.tms", &FileData); 
#endif

	if (hSearch == INVALID_HANDLE_VALUE) 
	{ 
	    return 0; // "No tms files found." 
	} 

	while (!fFinished) 
	{ 
#ifdef WIN_NT
		filenames[i] = new wchar_t[wcslen(FileData.cFileName)+1];
		wcscpy(filenames[i], FileData.cFileName);
#else 
		int iLength = strlen(FileData.cFileName);
		filenames[i] = new wchar_t[iLength];
		mbstowcs(filenames[i], FileData.cFileName, iLength);
		filenames[i][iLength] = 0;
#endif
	
#ifdef WIN_NT
		if (!FindNextFileW(hSearch, &FileData)) 
#else 
		if (!FindNextFileA(hSearch, &FileData)) 
#endif
		{
	        if (GetLastError() == ERROR_NO_MORE_FILES) 
			{ 
				fFinished = TRUE; 
			} 
			else 
			{ 
	            return i; //Couldn't find next file."; 
			} 
		}
		i++;
	}
	return i;
}

int TMOW32_API OpenLibrary(wchar_t* filename, TMO** pOperators)
{
	HMODULE hLibrary;
	map<wstring,HMODULE>::iterator i;
	int (*EnumOperators)(TMO**);
	wstring name(filename);

	i = mapModules.find(name);
	if (i == mapModules.end()) 
	{
#ifdef WIN_NT
		hLibrary = LoadLibraryW(filename);
#else 
		int iLength = wcslen(filename);
		char *filenameA = new char[iLength];
		wcstombs(filenameA, filename, iLength);
		filenameA[iLength] = 0;
		hLibrary = LoadLibraryA(filenameA);
#endif

		if (!hLibrary) return 0;
		mapModules.insert(map<wstring,HMODULE>::value_type(name, hLibrary));
	}
	else
	{
		hLibrary = (*i).second;
	}
	EnumOperators = (int(*)(TMO**)) GetProcAddress(hLibrary,"EnumOperators");
	if (EnumOperators) return EnumOperators(pOperators);
	return 0;
}

int TMOW32_API CloseLibrary(wchar_t* filename, TMO** pOperators)
{
	HMODULE hLibrary;
	map<wstring,HMODULE>::iterator i;
	int (*DeleteOperators)(TMO**);
	int retval = 0;

	wstring name(filename);

	i = mapModules.find(name);
	if (i == mapModules.end()) 
	{
		return -1;
	}
	else
	{
		hLibrary = (*i).second;
		DeleteOperators = (int(*)(TMO**)) GetProcAddress(hLibrary, "DeleteOperators");
		if (DeleteOperators) retval = DeleteOperators(pOperators);
		FreeLibrary(hLibrary);
		mapModules.erase(i);
	}
	return retval;	
}

int TMOW32_API OperatorCount(wchar_t* filename)
{
	HMODULE hLibrary;
	map<wstring,HMODULE>::iterator i;
	int (*OperatorCount)();

	wstring name(filename);

	i = mapModules.find(name);
	if (i == mapModules.end()) 
	{
#ifdef WIN_NT
		hLibrary = LoadLibraryW(filename);
#else 
		int iLength = wcslen(filename);
		char *filenameA = new char[iLength];
		wcstombs(filenameA, filename, iLength);
		filenameA[iLength] = 0;
		hLibrary = LoadLibraryA(filenameA);
#endif

		if (!hLibrary) return 0;
		mapModules.insert(map<wstring,HMODULE>::value_type(name, hLibrary));
	}
	else
	{
		hLibrary = (*i).second;
	}
	OperatorCount = (int(*)()) GetProcAddress(hLibrary, "OperatorCount");
	if (OperatorCount) return OperatorCount();
	return 0;
}