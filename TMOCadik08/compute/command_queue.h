#ifndef VOCCEL_COMPUTE_COMMAND_QUEUE_H
#define VOCCEL_COMPUTE_COMMAND_QUEUE_H

#include <vector>
#include "object.h"
#include "event.h"

namespace cl
{

class device;
class context;
class buffer;
class kernel;

typedef std::vector<size_t> ndrange;
typedef std::vector<event> event_list;

class command_queue : public object<command_queue, cl_command_queue> {
	public:
	command_queue(const device&, const context&,
	              const cl_command_queue_properties = 0);

	void retain() const;
	void release() const;
	template <bool = true>
	event read_buffer(const buffer&, const size_t, const size_t,
	                  void* const, const event_list = {}) const;
	template <bool = true>
	event write_buffer(const buffer&, const size_t, const size_t,
	                   const void* const, event_list = {}) const;
	event copy_buffer(const buffer&, const buffer&, const size_t,
	                  const size_t, const size_t, event_list = {}) const;
	template <bool = true>
	event map_buffer(void*, const buffer&, const cl_map_flags,
	                 const size_t, const size_t, const event_list = {}) const;
	//event unmap_buffer(void*, const buffer&, void* const,
	//                   const event_list = {}) const;
	event ndrange_kernel(const kernel&, const ndrange, ndrange,
	                     const ndrange, const event_list = {}) const;
	void wait(const event_list) const;
};

}// namespace cl

#include "compute/command_queue.inl"

#endif
