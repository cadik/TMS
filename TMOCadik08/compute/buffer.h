#ifndef VOCCEL_COMPUTE_BUFFER_H
#define VOCCEL_COMPUTE_BUFFER_H

#include "mem.h"

namespace cl
{

class context;

class buffer : public mem {
	public:
	buffer();
	buffer(const cl_mem&);
	buffer(const context&, const cl_mem_flags, const size_t,
	       const void* const = nullptr);

	template <typename T>
	cl::buffer create_sub(const cl_mem_flags, const size_t, const size_t,
	                      const unsigned, unsigned* const = nullptr) const;
};

}// namespace cl

#include "compute/buffer.inl"

#endif
