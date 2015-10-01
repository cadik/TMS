#include "TMOLinux.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

#include <sys/types.h>
#include <dirent.h>
#include <dlfcn.h>

#include <wchar.h>
#include <string>
#include "handle_map.h"

handle_map<std::string,void *> mapModules;

class TMO;

int EnumLibraries(wchar_t** filenames, int maxcount)
{
 DIR * d;
 struct dirent * de;
 d=opendir(".");
 int cnt=0;
 while((de=readdir(d)))
 {
  unsigned l = strlen(de->d_name);
  if(l>3 && !strncmp(de->d_name+(l-4),".tml",4))
  {
   if(cnt==maxcount)return cnt;
   filenames[cnt] = new wchar_t[l+1];
   mbstowcs(filenames[cnt], de->d_name, l);
   filenames[cnt][l] = 0;
   cnt++;
  }
 }

 return cnt;
}

void * get_function(wchar_t * filename,const char * function, void ** handle=NULL)
{
 unsigned l = wcslen(filename);
 char * fn = new char[l+3];
 wcstombs(fn+2,filename,l);
 fn[0]='.';
 fn[1]='/';
 fn[l+2]=0;
 std::string fns(fn+2);

 handle_map<std::string, void *>::iterator i = mapModules.find(fns);

 void * dlhandle=NULL;
 char * err;
 if(i==mapModules.end())
 {
  dlhandle = dlopen(fn,RTLD_LAZY);
  if(!dlhandle)
  {
   std::cerr << dlerror() << std::endl;
   return NULL;
  }
  mapModules[fns]=dlhandle;
 // std::cerr << "new dlhandle " << fns << " " << dlhandle << std::endl;
 }
 else
 {
  dlhandle=mapModules[fns];
 // std::cerr << "old dlhandle " << fns << " " << dlhandle << std::endl;
 }
 
 if(handle)*handle = dlhandle;

 void * rv = dlsym(dlhandle,function);
 if ((err = dlerror()) != NULL)
 {
  std::cerr << err << std::endl;
  dlclose(dlhandle);
  return NULL;
 }
 return rv;
}


int OpenLibrary(wchar_t* filename, TMO** pOperators)
{
 int (*EnumOperators)(TMO**);
 EnumOperators = (int(*)(TMO**))get_function(filename,"EnumOperators");
 if(EnumOperators)
 {
  return EnumOperators(pOperators);
 }
 else return 0;
}

int OperatorCount(wchar_t* filename)
{
 int (*OperatorCount)();
 OperatorCount = (int(*)())get_function(filename,"OperatorCount");
 if(OperatorCount)return OperatorCount();
 else return 0;
}

int CloseLibrary(wchar_t* filename, TMO** pOperators)
{
 void * dlhandle = NULL;
 int (*DeleteOperators)(TMO**);
 int rv = 0;
 DeleteOperators = (int(*)(TMO**))get_function(filename,"DeleteOperators",&dlhandle);
 if(DeleteOperators) rv = DeleteOperators(pOperators);
 if(dlhandle)dlclose(dlhandle);

 unsigned l = wcslen(filename);
 char * fn = new char[l+1];
 fn[l]=0;
 wcstombs(fn,filename,l);
 std::string fns(fn);
 mapModules.erase(fns);
 return rv;
}

