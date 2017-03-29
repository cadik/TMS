#ifndef VOCCEL_COMPUTE_PLATFORM_H
#define VOCCEL_COMPUTE_PLATFORM_H

#include "object.h"

namespace cl
{

class platform : public object<platform, cl_platform_id> {
	public:
	platform();

	void retain() const;
	void release() const;
};

}// namespace cl

#include "compute/platform.inl"

#endif
