#ifndef LQSTRING
#define LQSTRING

#include <qstring.h>

//QString setUnicodeCodes method doesn`t accept this string in wchar_t. This class provides an overloaded mathod to accept this and makes necessary conversion from ushort
class lqstring : public QString
{
 public:
#ifdef LINUX
  virtual QString & setUnicodeCodes(const wchar_t * s, size_t l);
#endif
  virtual ~lqstring(){}
};


#endif
