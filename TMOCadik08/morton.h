#ifndef TMOCADIK08_MORTON_H
#define TMOCADIK08_MORTON_H

typedef unsigned long long morton;

morton morton2d_encode(const unsigned x, const unsigned y);
void morton2d_decode(const morton m, unsigned& x, unsigned& y);

#endif
