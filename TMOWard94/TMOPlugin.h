#include "TMO.h"

#ifndef LINUX
 #ifdef TMOPLUGIN_EXPORTS
 #define TMOPLUGIN_API __declspec(dllexport)
 #else
 #define TMOPLUGIN_API __declspec(dllimport)
 #endif
 #pragma warning (disable: 4251) 
#else
 #define TMOPLUGIN_API
#endif


extern "C" TMOPLUGIN_API int EnumOperators(TMO** operators);
extern "C" TMOPLUGIN_API int DeleteOperators(TMO** operators);
extern "C" TMOPLUGIN_API int OperatorCount();
