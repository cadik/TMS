#ifndef COMPUTE_PLATFORM_H
#define COMPUTE_PLATFORM_H

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

#include "platform.inl"

#endif
