#include "common/error.h"

namespace com
{

inline loader::loader(const std::string& path) : input{path}
{
	if (!input.is_open())
		throw error{__FILE__, __LINE__, __PRETTY_FUNCTION__, "std::ifstream::is_open: input error on \"" + path + '\"'};
}

inline loader::~loader()
{
	input.close();
}

}// namespace com
