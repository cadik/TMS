#include "lqstring.h"

#ifdef LINUX

//convert string in wchar_t array to ushort array
ushort * wchar2ushort(const wchar_t * s)
{
 unsigned l = wcslen(s);
 ushort * us = new ushort[l+1];
 ushort * rv = us;
 us[l]=0;
 while(*s)
 {
  *us = *s;
   us++;
   s++;
 }
 return rv;
}


//create lqstring containing the string pointed by s
lqstring GetString(const wchar_t* s)
{
  lqstring str;
  ushort * us = wchar2ushort(s);
  str.QString::setUnicodeCodes(us, wcslen(s));
  delete[] us;
  return str;
}



//the overloaded method
QString & lqstring::setUnicodeCodes(const wchar_t * s, size_t l)
{
 (*this)=GetString(s);
 return (*this);
}

#endif

