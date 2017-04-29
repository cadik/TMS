#ifndef VOCCEL_COMPUTE_PARALLEL_H
#define VOCCEL_COMPUTE_PARALLEL_H

#include "compute/platform.h"
#include "compute/device.h"
#include "compute/context.h"
#include "compute/command_queue.h"
#include "compute/buffer.h"
#include "compute/program.h"

namespace cl
{

class parallel {
	public:
	parallel(const cl_device_type = CL_DEVICE_TYPE_CPU);

	command_queue create_command_queue(const cl_command_queue_properties = 0) const;
	buffer create_buffer(const cl_mem_flags, const size_t,
	                     const void* const = nullptr) const;
	program create_program(const char*, const char* const = nullptr) const;
	size_t get_wgs() const;
	unsigned get_maba() const;
	std::string info() const;

	private:
	const cl::platform host;
	const cl::device gpu;
	const cl::context env;
	const size_t lms, // local mem. size
	             wgs; // work group size
	const unsigned maba; // memory addr. base align
};

}// namespace cl

#include "compute/parallel.inl"

#endif
