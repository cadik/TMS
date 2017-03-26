#ifndef VOCCEL_COMPUTE_EVENT_H
#define VOCCEL_COMPUTE_EVENT_H

#include "object.h"

namespace cl
{

class event : public object<event, cl_event> {
	public:
	void retain() const;
	void release() const;

	void profiling_info(const cl_profiling_info, const size_t,
	                    void* const) const;
};

}// namespace cl

#include "compute/event.inl"

#endif
