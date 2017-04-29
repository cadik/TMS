#ifndef TMOCADIK08_SPACEFILL_H
#define TMOCADIK08_SPACEFILL_H

#include <vector>
#include "vec.h"

std::vector<vec2i> hilbert(const unsigned n, const float x0, const float y0,
                           const float xi, const float xj, const float yi,
                           const float yj);
std::vector<vec2i> moore(const unsigned n, const float x0, const float y0,
                         const float xi, const float xj, const float yi,
                         const float yj);

#endif
