#include <iostream>
#include "common/error.h"

namespace cl
{

inline context::context(const properties opts, const device& dev)
{
	const auto pfn_notify{[] (const char* errinfo, const void*, size_t, void*) {
		std::cerr << errinfo << std::endl;
	}};
	cl_int ret;
	id = clCreateContext(opts.begin(), 1, &((const cl_device_id&) dev), pfn_notify, nullptr, &ret);
	if (ret != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "clCreateContext: " + std::to_string(ret)};
}

inline void context::retain() const
{
	clRetainContext(id);
}

inline void context::release() const
{
	clReleaseContext(id);
}

}// namespace cl
