#ifndef COMMON_LOADER_H
#define COMMON_LOADER_H

#include <string>
#include <fstream>

namespace com
{

struct loader {
	loader(const std::string&);
	~loader();

	protected:
	std::ifstream input;
};

}// namespace com

#include "loader.inl"

#endif
