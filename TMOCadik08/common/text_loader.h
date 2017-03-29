#ifndef COMMON_TEXT_LOADER_H
#define COMMON_TEXT_LOADER_H

#include <streambuf>
#include "loader.h"

namespace com
{

struct text_loader : public loader {
	text_loader(const std::string);

	std::string operator()();
};

}// namespace com

#include "text_loader.inl"

#endif
