//
// Created by daeyun on 5/14/17.
//

#ifndef MESHDIST_MESHDIST_CGAL_H
#define MESHDIST_MESHDIST_CGAL_H

#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <omp.h>

#include "meshdist.h"

namespace meshdist_cgal {
using meshdist::Triangle;

float MeshToMeshDistanceOneDirection(const std::vector<Triangle> &from,
                                     const std::vector<Triangle> &to,
                                     float sampling_density);

float MeshToMeshDistance(const std::vector<Triangle> &a, const std::vector<Triangle> &b);
}

#endif //MESHDIST_MESHDIST_CGAL_H
