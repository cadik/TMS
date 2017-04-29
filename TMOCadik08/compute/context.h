#ifndef VOCCEL_COMPUTE_CONTEXT_H
#define VOCCEL_COMPUTE_CONTEXT_H

#include "object.h"

namespace cl
{

class platform;
class device;

class context : public object<context, cl_context> {
	typedef std::initializer_list<cl_context_properties> properties;

	public:
	context(const properties, const device&);

	void retain() const;
	void release() const;
};

}// namespace cl

#include "compute/context.inl"

#endif
