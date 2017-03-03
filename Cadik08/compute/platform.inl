#include "common/error.h"

namespace cl
{

inline platform::platform()
{
	if (const cl_int ret = clGetPlatformIDs(1, &id, nullptr) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "clGetPlatformIDs: " + std::to_string(ret)};
}

inline void platform::retain() const
{
} 

inline void platform::release() const
{
}

}// namespace cl
