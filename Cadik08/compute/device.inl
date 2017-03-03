#include "common/error.h"

namespace cl
{

//______________________________________________________________________________
inline device::device(const platform& host, const cl_device_type type)
{
	if (const cl_int ret = clGetDeviceIDs((const cl_platform_id&) host,
	                                      type, 1, &id, nullptr) !=
	                                      CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clGetDeviceIDs: " + std::to_string(ret)};
}

//______________________________________________________________________________
inline void device::retain() const
{
} 

//______________________________________________________________________________
inline void device::release() const
{
}

//______________________________________________________________________________
template <typename T>
inline T device::info(const cl_device_info name) const
{
	//size_t size;
	//clGetDeviceInfo(id, name, 0, nullptr, &size);
	T val;
	if (const cl_int ret = clGetDeviceInfo(id, name, sizeof(T), &val, nullptr))
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clGetDeviceInfo: " + std::to_string(ret)};
	return val;
}

//______________________________________________________________________________
template <>
inline char* device::info(const cl_device_info name) const
{
	size_t size;
	clGetDeviceInfo(id, name, 0, nullptr, &size);
	char* vals{new char[size]};
	if (const cl_int ret = clGetDeviceInfo(id, name, size, vals, nullptr))
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clGetDeviceInfo: " + std::to_string(ret)};
	return vals;
}

}// namespace cl
