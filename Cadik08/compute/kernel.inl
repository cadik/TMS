#include "common/error.h"

namespace cl
{

class buffer;

inline kernel::kernel()
{
}

inline kernel::kernel(const program& exe, const char* const name)
{
	cl_int ret;
	id = clCreateKernel((const cl_program&) exe, name, &ret);
	if (ret != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "clCreateKernel: " + std::to_string(ret)};
}

inline void kernel::retain() const
{
	clRetainKernel(id);
}

inline void kernel::release() const
{
	clReleaseKernel(id);
}

template <typename T>
inline void kernel::set_arg(const cl_uint i, const T& val) const
{
	if (const cl_int ret = clSetKernelArg(id, i, sizeof(T), &val) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "clSetKernelArg(" + std::to_string(i) + "): " + std::to_string(ret)};
}

template <>
inline void kernel::set_arg<buffer>(const cl_uint i, const buffer& val) const
{
	if (const cl_int ret = clSetKernelArg(id, i, sizeof(cl_mem), (const cl_mem*) val) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "clSetKernelArg(" + std::to_string(i) + "): " + std::to_string(ret)};
}

template <>
inline void kernel::set_arg<local_mem>(const cl_uint i, const local_mem& val) const
{
	if (const cl_int ret = clSetKernelArg(id, i, val.n, nullptr) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "clSetKernelArg(" + std::to_string(i) + "): " + std::to_string(ret)};
}

template <typename ...Ts>
inline void kernel::set_args(const Ts&... args) const
{
	cl_uint i{};
	std::initializer_list<int>{(set_arg(i++, args), 0)...};
}

}// namespace cl
