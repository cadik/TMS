#ifndef VOCCEL_COMPUTE_PLATFORM_H
#define VOCCEL_COMPUTE_PLATFORM_H

#include "object.h"

namespace cl
{

class platform : public object<platform, cl_platform_id> {
	public:
	platform(const cl_device_type);

	void retain() const;
	void release() const;
};

}// namespace cl

#include "compute/platform.inl"

#endif
