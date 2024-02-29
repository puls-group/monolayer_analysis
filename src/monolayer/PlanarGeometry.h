#ifndef __PERIODIC_GEOMETRY_H__
#define __PERIODIC_GEOMETRY_H__
#include <cmath>
#include <iostream>

#include <complex>

using v2d = std::complex<double>;

// Restrict vector to parallel epipede spanned by pbc. 
// Maximum coefficient: +-0.5 in any direction after restriction
v2d restrict(const v2d vector, const std::pair<v2d, v2d> &pbc);

// Represent vector in terms of basis given by base
v2d represent(const v2d vector, const std::pair<v2d, v2d> &base);

// Does a1->b1 intersect segment a2->b2
bool intersectsSegment(v2d a1, v2d b1, v2d a2, v2d b2);

#endif
