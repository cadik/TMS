#include <cmath>
#include "quadtree.h"
#include "common/math.h"

//______________________________________________________________________________
quadtree::quadtree(const unsigned n) : side{com::math::ceil2pow(n, 4)},
                                       height{std::log(side * side) /
                                              std::log(4)},
                                       root((4 * side * side - 1) / 3),
                                       offset0{root.size() - std::pow(4, height)}
{
}

//______________________________________________________________________________
const vec2d& quadtree::operator[](const unsigned i) const
{
	return root[offset0 + i];
} 

//______________________________________________________________________________
vec2d& quadtree::operator[](const unsigned i)
{
	return root[offset0 + i];
} 

//______________________________________________________________________________
unsigned quadtree::get_side() const
{
	return side;
}

//______________________________________________________________________________
unsigned quadtree::get_height() const
{
	return height;
}

//______________________________________________________________________________
unsigned quadtree::size() const
{
	return root.size();
}

//______________________________________________________________________________
unsigned quadtree::len() const
{
	return root.size() - offset0;
}

//______________________________________________________________________________
vec2d* quadtree::data()
{
	return root.data();
}

//______________________________________________________________________________
std::string quadtree::print() const
{
	std::string info{};

	for (unsigned i = 0; i < height; ++i) {
		const unsigned offset_i = (4 * pow(4, i) - 1) / 3,
			       size_i = ((4 * pow(4, i + 1) - 1) / 3) -
		                        offset_i;
		info += "-----------------\n";
		info += std::to_string(offset_i) + ' ' + std::to_string(size_i) + '\n';

		/*if (i < 3) {
		for (unsigned i = offset_i; i < size_i + offset_i; ++i)
			info += '(' + std::to_string(root[i].x) + ", " + std::to_string(root[i].y) + "), ";
		info += '\n';
		}*/
	}

	return info;
}
