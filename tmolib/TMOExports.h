#ifndef LINUX

 #ifdef TMOLIB_EXPORTS
 #define TMOLIB_API __declspec(dllexport)
 #else
 #define TMOLIB_API __declspec(dllimport)
 #endif
 #pragma warning (disable: 4251) 

#else

 #define TMOLIB_API

#endif
