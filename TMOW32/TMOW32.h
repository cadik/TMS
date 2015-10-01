#ifdef TMOW32_EXPORTS
#define TMOW32_API __declspec(dllexport)
#else
#define TMOW32_API __declspec(dllimport)
#endif
#pragma warning (disable: 4251) 
#pragma warning (disable: 4786)
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

#include <windows.h>
class TMO;

int TMOW32_API CloseLibrary(wchar_t* filename, TMO** pOperators);
int TMOW32_API OpenLibrary(wchar_t* filename, TMO** pOperators);
int TMOW32_API EnumLibraries(wchar_t** filenames, int maxcount);
int TMOW32_API OperatorCount(wchar_t* filename);