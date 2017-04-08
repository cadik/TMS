#ifndef TMOCADIK08_QUADTREE_H
#define TMOCADIK08_QUADTREE_H

#include <vector>
#include <string>
#include "vec.h"

class quadtree {
	public:
	quadtree(const unsigned n);

	const vec2d& operator[](const unsigned) const;
	vec2d& operator[](const unsigned);

	unsigned get_side() const;
	unsigned get_height() const;
	unsigned size() const;
	unsigned len() const;
	vec2d* data();
	std::string print() const;

	private:
	const unsigned side,
	               height; // root node is not accounted for!
	std::vector<vec2d> root;
	const unsigned offset0;
};

#endif
