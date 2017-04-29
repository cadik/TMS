namespace com
{
namespace math
{

//______________________________________________________________________________
inline unsigned ceil2mul(const unsigned x, const unsigned short m)
{
	if (!m)
		return x;

	const unsigned short r = x % m;
	return r ? x + m - r : x;
}

//______________________________________________________________________________
inline unsigned floor2mul(const unsigned x, const unsigned short m)
{
	return m ? x - x % m : x;
}

//______________________________________________________________________________
inline unsigned ceil2pow(const unsigned x, const unsigned short e)
{
	unsigned pow{1};
	while (pow < x)
		pow *= e;
	return pow;
}

//______________________________________________________________________________
inline unsigned ceil2pow2(unsigned x)
{
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

}// namespace math
}// namespace com
