#include "spacefill.h"

//______________________________________________________________________________
std::vector<vec2i> hilbert(const unsigned n, const float x0, const float y0,
                           const float xi, const float xj, const float yi,
                           const float yj)
{
	std::vector<vec2i> pts{};
	if (!n)
		pts.push_back(vec2i{(int) (x0 + (xi + yi) / 2.f),
		                    (int) (y0 + (xj + yj) / 2.f)});
	else {
		const std::vector<vec2i> a = hilbert(n - 1,
		                                     x0, y0,
		                                     yi / 2.f, yj / 2.f,
		                                     xi / 2.f, xj / 2.f);
		const std::vector<vec2i> b = hilbert(n - 1,
		                                     x0 + xi / 2.f, y0 + xj / 2.f,
		                                     xi / 2.f, xj / 2.f,
		                                     yi / 2.f, yj / 2.f);
		const std::vector<vec2i> c = hilbert(n - 1,
		                                     x0 + xi / 2.f + yi / 2.f, y0 + xj / 2.f + yj / 2.f,
		                                     xi / 2.f, xj / 2.f,
		                                     yi / 2.f, yj / 2.f);
		const std::vector<vec2i> d = hilbert(n - 1,
		                                     x0 + xi / 2.f + yi, y0 + xj / 2.f + yj,
		                                     -yi / 2.f, -yj / 2.f,
		                                     -xi / 2.f, -xj / 2.f);
		std::copy(a.begin(), a.end(), std::back_inserter(pts));
		std::copy(b.begin(), b.end(), std::back_inserter(pts));
		std::copy(c.begin(), c.end(), std::back_inserter(pts));
		std::copy(d.begin(), d.end(), std::back_inserter(pts));
	}

	return pts;
}

//______________________________________________________________________________
std::vector<vec2i> moore(const unsigned n, const float x0, const float y0,
                         const float xi, const float xj, const float yi,
                         const float yj)
{
	std::vector<vec2i> pts{};
	if (!n)
		pts.push_back(vec2i{(int) (x0 + (xi + yi) / 2.f),
		                    (int) (y0 + (xj + yj) / 2.f)});
	else {
		const std::vector<vec2i> a = hilbert(n - 1,
		                                     x0 + xi / 2.f, y0 + xj / 2.f,
		                                     -xi / 2.f, xj / 2.f,
		                                     yi / 2.f, yj / 2.f);
		const std::vector<vec2i> b = hilbert(n - 1,
		                                     x0 + xi / 2.f + yi / 2.f, y0 + xj / 2.f + yj / 2.f,
		                                     -xi / 2.f, xj / 2.f,
		                                     yi / 2.f, yj / 2.f);
		const std::vector<vec2i> c = hilbert(n - 1,
		                                     x0 + xi / 2.f + yi, y0 + xj / 2.f + yj,
		                                     xi / 2.f, xj / 2.f,
		                                     yi / 2.f, -yj / 2.f);
		const std::vector<vec2i> d = hilbert(n - 1,
		                                     x0 + xi / 2.f + yi / 2.f, y0 + xj / 2.f + yj / 2.f,
		                                     xi / 2.f, xj / 2.f,
		                                     yi / 2.f, -yj / 2.f);
		std::copy(a.begin(), a.end(), std::back_inserter(pts));
		std::copy(b.begin(), b.end(), std::back_inserter(pts));
		std::copy(c.begin(), c.end(), std::back_inserter(pts));
		std::copy(d.begin(), d.end(), std::back_inserter(pts));
	}

	return pts;
}
