#include "TMO.h"

// Just a drop-in replacement that ChatGPT suggested, otherwise it wouldn't compile
#ifdef _WIN32
#define TMOPLUGIN_API __declspec(dllexport)
#else
#define TMOPLUGIN_API __attribute__((visibility("default")))
#endif

extern "C" TMOPLUGIN_API int EnumOperators(TMO **operators);
extern "C" TMOPLUGIN_API int DeleteOperators(TMO **operators);
extern "C" TMOPLUGIN_API int OperatorCount();