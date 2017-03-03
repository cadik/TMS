#ifndef COMPUTE_PROGRAM_H
#define COMPUTE_PROGRAM_H

#include <map>
#include <string>
#include "object.h"
#include "kernel.h"

namespace cl
{

class device;
class context;

class program : public object<program, cl_program> {
	public:
	program(const context&, const char*, const device&, const char* const = nullptr);

	kernel& operator[](const std::string);

	void retain() const;
	void release() const;

	private:
	std::map<std::string, kernel> tasks;
};

}// namespace cl

#include "program.inl"

#endif
