#ifndef __GRAPH_CONSTRUCTOR_H__
#define __GRAPH_CONSTRUCTOR_H__
#include "Graph2d.h"
#include <vector>

using v2d = std::complex<double>;
Graph2D *constructPlanarGraph(const std::vector<v2d> &points, const std::pair<v2d, v2d> &pbc, double cutoff_dist);

#endif
