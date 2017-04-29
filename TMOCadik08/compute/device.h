#ifndef VOCCEL_COMPUTE_DEVICE_H
#define VOCCEL_COMPUTE_DEVICE_H

#include "object.h"

namespace cl
{

class platform;

class device : public object<device, cl_device_id> {
	public:
	device(const platform&, const cl_device_type);

	void retain() const;
	void release() const;

	template <typename T>
	T info(const cl_device_info) const;
};

}// namespace cl

#include "compute/device.inl"

#endif
