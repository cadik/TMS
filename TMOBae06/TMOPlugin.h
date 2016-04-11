#include "TMO.h"
/*
#ifdef TMOPLUGIN_EXPORTS
#define TMOPLUGIN_API __declspec(dllexport)
#else
#define TMOPLUGIN_API __declspec(dllimport)
#endif
#pragma warning (disable: 4251) 
*/

#define TMOPLUGIN_API

extern "C" TMOPLUGIN_API int EnumOperators(TMO** operators);
extern "C" TMOPLUGIN_API int DeleteOperators(TMO** operators);
extern "C" TMOPLUGIN_API int OperatorCount();

