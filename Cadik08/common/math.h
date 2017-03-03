#ifndef COMMON_MATH_H
#define COMMON_MATH_H

namespace com
{
namespace math
{

unsigned ceil_multiple(const unsigned, const unsigned short);
unsigned floor_multiple(const unsigned, const unsigned short);
unsigned ceil_power(const unsigned, const unsigned short);

}// namespace math
}// namespace com

#include "math.inl"

#endif
