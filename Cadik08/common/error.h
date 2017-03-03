#ifndef COMMON_ERROR_H
#define COMMON_ERROR_H

#include <stdexcept>
#include <string>

namespace com
{

class error : public std::runtime_error {
	public:
	error(const std::string, const unsigned, const char* const, const std::string);
};

}// namespace com

#include "error.inl"

#endif
