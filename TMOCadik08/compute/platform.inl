#include "common/error.h"

namespace cl
{

//______________________________________________________________________________
inline platform::platform()
{
	/*cl_platform_id ids[2];
	clGetPlatformIDs(2, ids, nullptr);
	id = ids[0];*/
	if (const cl_int ret = clGetPlatformIDs(1, &id, nullptr) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clGetPlatformIDs: " + std::to_string(ret)};
}

//______________________________________________________________________________
inline void platform::retain() const
{
} 

//______________________________________________________________________________
inline void platform::release() const
{
}

}// namespace cl
