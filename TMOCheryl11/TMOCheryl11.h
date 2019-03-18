#include "TMO.h"

class TMOCheryl11 : public TMO
{
public:
    TMOCheryl11();
    virtual ~TMOCheryl11();
    virtual int Transform();

protected:
    TMODouble dParameter;

    void clusterize(bool showClusteredImg);
};
