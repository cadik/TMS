#include "common/error.h"
#include "common/math.h"

namespace cl
{

//______________________________________________________________________________
inline command_queue::command_queue(const device& dev, const context& env,
                                    const cl_command_queue_properties opts)
{
	cl_int ret;
	id = clCreateCommandQueue((const cl_context&) env,
	                          (const cl_device_id&) dev, opts, &ret);
	if (ret != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clCreateCommandQueue: " +
		                 std::to_string(ret)};
}

//______________________________________________________________________________
inline void command_queue::retain() const
{
	clRetainCommandQueue(id);
}

//______________________________________________________________________________
inline void command_queue::release() const
{
	clFlush(id);
	clFinish(id);
	clReleaseCommandQueue(id);
}

//______________________________________________________________________________
template <bool B>
inline event command_queue::read_buffer(const buffer& buf, const size_t offset,
                                        const size_t cb, void* const ptr,
                                        const event_list pending) const
{
	event status;
	if (const cl_int ret = clEnqueueReadBuffer(id, (const cl_mem&) buf, B,
	                                           offset, cb, ptr,
	                                           pending.size(),
	                                           (const cl_event*) pending[0],
	                                           (cl_event*) status) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clEnqueueReadBuffer: " + std::to_string(ret)};
	return status;
}

//______________________________________________________________________________
template <bool B>
inline event command_queue::write_buffer(const buffer& buf, const size_t offset,
                                         const size_t cb, const void* const ptr,
                                         const event_list pending) const
{
	event status;
	if (const cl_int ret = clEnqueueWriteBuffer(id, (const cl_mem&) buf, B,
	                                            offset, cb, ptr,
	                                            pending.size(),
	                                            (const cl_event*) pending[0],
	                                            (cl_event*) status) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clEnqueueWriteBuffer: " + std::to_string(ret)};
	return status;
}

//______________________________________________________________________________
inline event command_queue::copy_buffer(const buffer& src, const buffer& dst,
                                        const size_t src_ofs,
                                        const size_t dst_ofs, const size_t cb,
                                        event_list pending) const
{
	event status;
	if (const cl_int ret = clEnqueueCopyBuffer(id, (const cl_mem&) src,
	                                           (const cl_mem&) dst, src_ofs,
	                                           dst_ofs, cb, pending.size(),
	                                           (const cl_event*) pending[0],
	                                           (cl_event*) status) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clEnqueueWriteBuffer: " + std::to_string(ret)};
	return status;
}

//______________________________________________________________________________
template <bool B>
inline event command_queue::map_buffer(void* ptr, const buffer& buf,
                                       const cl_map_flags opts,
                                       const size_t offset, const size_t cb,
                                       const event_list pending) const
{
	event status;

	cl_int ret;
	ptr = clEnqueueMapBuffer(id, (const cl_mem&) buf, B, opts,
	                         offset, cb, pending.size(),
	                         (const cl_event*) &pending[0], (cl_event*) status,
	                         &ret);
	if (ret != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clEnqueueMapBuffer: " + std::to_string(ret)};

	return status;
}

//______________________________________________________________________________
inline event command_queue::ndrange_kernel(const kernel& task,
                                           const ndrange offset,
                                           ndrange global,
                                           const ndrange local,
                                           const event_list pending) const
{
	for (size_t i{}; i < local.size(); ++i)
		global[i] = com::math::ceil2mul(global[i], local[i]);

	event status;
	if (const cl_int ret = clEnqueueNDRangeKernel(id,
	                                              (const cl_kernel&) task,
	                                              global.size(), &offset[0],
	                                              &global[0], &local[0],
	                                              pending.size(),
	                                              (const cl_event*) pending[0],
	                                              (cl_event*) status) != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__,
		                 "clEnqueueNDRangeKernel: " +
		                 std::to_string(ret)};
	return status;
}

//______________________________________________________________________________
inline void command_queue::wait(const event_list pending) const
{
	clWaitForEvents(pending.size(), (const cl_event*) pending.data());
}

}// namespace cl
