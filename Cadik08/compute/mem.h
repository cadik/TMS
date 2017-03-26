#ifndef VOCCEL_COMPUTE_MEM_H
#define VOCCEL_COMPUTE_MEM_H

#include "object.h"

namespace cl
{

class mem : public object<mem, cl_mem> {
	public:
	mem();
	mem(const cl_mem&);

	void retain() const;
	void release() const;
};

}// namespace cl

#include "compute/mem.inl"

#endif
