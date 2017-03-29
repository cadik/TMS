#include "common/error.h"

namespace cl
{

//______________________________________________________________________________
inline void event::retain() const
{
	clRetainEvent(id);
}

//______________________________________________________________________________
inline void event::release() const
{
	clReleaseEvent(id);
}

//______________________________________________________________________________
inline void event::profiling_info(const cl_profiling_info param,
                                  const size_t size, void* const val) const
{
	if (const cl_int ret = clGetEventProfilingInfo(id, param, size, val,
	                                               nullptr) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clGetEventProfilingInfo: " + std::to_string(ret)};
}

}// namespace cl
