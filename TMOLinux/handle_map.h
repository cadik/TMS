#include <map>


template <class S,class T> class handle_map: public std::map<S,T>
{
 public:
 typedef typename std::map<S,T>::iterator it;
  ~handle_map()
  {
   it i = this->begin();
   while(i!=this->end())
   {
    dlclose(i->second);
    i++;
   }
  }
};
