#ifndef VOCCEL_COMPUTE_KERNEL_H
#define VOCCEL_COMPUTE_KERNEL_H

#include "object.h"

namespace cl
{

class program;

struct local_mem {
	size_t n;
};

class kernel : public object<kernel, cl_kernel> {
	public:
	kernel();
	kernel(const program&, const char* const);

	void retain() const;
	void release() const;

	template<typename T>
	void set_arg(const cl_uint, const T&) const;
	template<typename ...Ts>
	void set_args(const Ts&...) const;
};

}// namespace cl

#include "compute/kernel.inl"

#endif
