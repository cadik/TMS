#include "TMO.h"
#include <torch/torch.h>
#include <torch/script.h>
#include <opencv2/opencv.hpp>

class TMOJah06 : public TMO
{
public:
	TMOJah06();
	virtual ~TMOJah06();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
