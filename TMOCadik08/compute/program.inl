#include <string>
#include "common/error.h"

namespace cl
{

inline program::program(const context& env, const char* source, const device& dev, const char* const opts)
{
	cl_int ret;
	id = clCreateProgramWithSource((const cl_context&) env, 1, &source, nullptr, &ret);
	if (ret != CL_SUCCESS)
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "clCreateProgramWithSource: " + std::to_string(ret)};

	if ((ret = clBuildProgram(id, 1, &((const cl_device_id&) dev), opts, nullptr, nullptr) != CL_SUCCESS)) {
		std::string msg; 
		//if (ret == CL_BUILD_PROGRAM_FAILURE) {
			size_t log_size;
			clGetProgramBuildInfo(id, (const cl_device_id&) dev, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
			char* log{new char[log_size]};
			clGetProgramBuildInfo(id, (const cl_device_id&) dev, CL_PROGRAM_BUILD_LOG, log_size, log, nullptr);
			msg += ": compiler log:\n";
			msg += log;
			delete[] log;
		//}
		throw com::error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "clBuildProgram: " + std::to_string(ret) + msg};
	}
}

inline kernel& program::operator[](const std::string name)
{
	if (tasks.find(name) == tasks.end())
		tasks.insert({name, {*this, name.c_str()}});
	return tasks[name];
}

inline void program::retain() const
{
	clRetainProgram(id);
}

inline void program::release() const
{
	clReleaseProgram(id);
}

}// namespace cl
