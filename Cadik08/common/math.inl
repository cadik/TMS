namespace com
{
namespace math
{

//______________________________________________________________________________
inline unsigned ceil_multiple(const unsigned x, const unsigned short m)
{
	if (!m)
		return x;

	const unsigned short r = x % m;
	return r ? x + m - r : x;
}

//______________________________________________________________________________
inline unsigned floor_multiple(const unsigned x, const unsigned short m)
{
	return m ? x - x % m : x;
}

//______________________________________________________________________________
inline unsigned ceil_power(const unsigned x, const unsigned short e)
{
	unsigned pow{1};
	while (pow < x)
		pow *= e;
	return pow;
}

}// namespace math
}// namespace com
