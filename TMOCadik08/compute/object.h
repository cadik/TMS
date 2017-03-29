#ifndef VOCCEL_COMPUTE_OBJECT_H
#define VOCCEL_COMPUTE_OBJECT_H

#include <CL/cl.h>

namespace cl
{

template <class D, typename T>
class object {
	public:
	object();
	object(const object&);
	object(const T&);
	~object();

	operator const T&() const;
	operator T&();
	operator const T*() const;
	operator T*();
	object& operator=(const object&);

	protected:
	T id;

	private:
	void retain() const;
	void release() const;
};

}// namespace cl

#include "compute/object.inl"

#endif
