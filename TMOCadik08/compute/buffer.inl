#include "common/error.h"
#include "common/math.h"

namespace cl
{

//______________________________________________________________________________
inline buffer::buffer()
{
}

//______________________________________________________________________________
inline buffer::buffer(const cl_mem& id) : mem{id}
{
}

//______________________________________________________________________________
inline buffer::buffer(const context& dev_ctx, const cl_mem_flags opts,
                      const size_t size, const void* const ptr)
{
	cl_int ret;
	id = clCreateBuffer((const cl_context&) dev_ctx, opts, size,
	                    const_cast<void* const>(ptr), &ret);
	if (ret != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clCreateBuffer: " + std::to_string(ret)};
}

//______________________________________________________________________________
template <typename T>
inline cl::buffer buffer::create_sub(const cl_mem_flags opts, const size_t at,
                                     const size_t amount, const unsigned maba,
                                     unsigned* const align) const
{
	const size_t origin = com::math::floor2mul(at * sizeof(T), maba),
		     offset = (at * sizeof(T) - origin) / sizeof(T);
	if (align)
		*align = offset;

	cl_int ret;
	const cl_buffer_region info{origin, (amount + offset) * sizeof(T)};
	const buffer sub{clCreateSubBuffer((const cl_mem&) id, opts,
	                 CL_BUFFER_CREATE_TYPE_REGION, &info, &ret)};
	if (ret != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clCreateSubBuffer: " + std::to_string(ret)};
	return sub;
}

}// namespace cl
