#include "common/error.h"

namespace cl
{

//______________________________________________________________________________
inline platform::platform(const cl_device_type type)
{
	cl_uint num_hosts = 0;
	clGetPlatformIDs(0, nullptr, &num_hosts);
	cl_platform_id* ids = new cl_platform_id[num_hosts];
	if (const cl_int ret = clGetPlatformIDs(num_hosts, ids, nullptr) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clGetPlatformIDs: " + std::to_string(ret)};

	for (unsigned i = 0; i < num_hosts; ++i) {
		cl_uint num_gpus = 0;
		clGetDeviceIDs((const cl_platform_id&) ids[i],
	                       type, 0, nullptr, &num_gpus); 
		if (num_gpus == 0)
			continue;

		id = ids[i];
		break;

	}
	delete[] ids;
	/*if (const cl_int ret = clGetPlatformIDs(1, &id, nullptr) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clGetPlatformIDs: " + std::to_string(ret)};*/
	/*cl_platform_id ids[2];
	clGetPlatformIDs(2, ids, nullptr);
	id = ids[1];*/
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
