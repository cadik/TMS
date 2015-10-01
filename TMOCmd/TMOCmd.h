// TMOCmd.h: interface for the TMOCmd class.
//
//////////////////////////////////////////////////////////////////////
#include "../tmolib/TMO.h"

class TMOCmd  
{
public:
	virtual int main(int argc, char *argv[]);
	TMOCmd();
	virtual ~TMOCmd();
protected:
	void Help(TMO** operators, int opindex = -1);
	static int Bar(int part, int all);

	int iOperatorCount;
};

