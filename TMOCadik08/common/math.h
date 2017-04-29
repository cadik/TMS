#ifndef COMMON_MATH_H
#define COMMON_MATH_H

namespace com
{
namespace math
{

unsigned ceil2mul(const unsigned, const unsigned short);
unsigned floor2mul(const unsigned, const unsigned short);
unsigned ceil2pow(const unsigned, const unsigned short);
unsigned ceil2pow2(unsigned);

}// namespace math
}// namespace com

#include "math.inl"

#endif
