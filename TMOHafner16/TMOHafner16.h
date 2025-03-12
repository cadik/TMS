#include "TMO.h"

class TMOHafner16 : public TMO
{
public:
TMOHafner16();
	virtual ~TMOHafner16();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
