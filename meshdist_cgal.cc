//
// Created by daeyun on 5/14/17.
//

#include "meshdist_cgal.h"

namespace meshdist_cgal {

typedef CGAL::Simple_cartesian<float> K;

typedef std::list<K::Triangle_3>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

using meshdist::Triangle;
using meshdist::Points3;

float MeshToMeshDistanceOneDirection(const std::vector<Triangle> &from,
                                     const std::vector<Triangle> &to,
                                     float sampling_density) {

  Points3 points;

  SamplePointsOnTriangles(from, sampling_density, &points);

  std::list<K::Triangle_3> triangle_list;
  for (const auto &triangle : to) {
    triangle_list.emplace_back(K::Point_3{triangle.a[0], triangle.a[1], triangle.a[2]},
                               K::Point_3{triangle.b[0], triangle.b[1], triangle.b[2]},
                               K::Point_3{triangle.c[0], triangle.c[1], triangle.c[2]});
  }
  std::vector<K::Point_3> point_list;
  for (int i = 0; i < points.cols(); ++i) {
    point_list.emplace_back(points(0, i), points(1, i), points(2, i));
  }

  int num_triangles = static_cast<int>(to.size());
  int num_points = static_cast<int>(points.cols());

  BOOST_LOG_TRIVIAL(debug) << "Computing minimum distances from " << num_points
                           << " points to " << num_triangles << " triangles.";

  auto start = meshdist::MilliSecondsSinceEpoch();
  Tree tree(triangle_list.begin(), triangle_list.end());
  tree.build();
  tree.accelerate_distance_queries();
  BOOST_LOG_TRIVIAL(debug) << "Time elapsed for building tree (CGAL): " << meshdist::MilliSecondsSinceEpoch() - start;

  float distance_sum = 0;

#pragma omp parallel for if(USE_OMP) reduction(+:distance_sum) schedule(SCHEDULE)
  for (int i = 0; i < point_list.size(); ++i) {
    float dist = tree.squared_distance(point_list[i]);
    distance_sum += dist;
  }

  BOOST_LOG_TRIVIAL(debug) << "distance: " << distance_sum;
  float rms = static_cast<float>(std::sqrt(distance_sum / static_cast<double>(point_list.size())));
  BOOST_LOG_TRIVIAL(debug) << "RMS: " << rms;
  auto elapsed = meshdist::MilliSecondsSinceEpoch() - start;
  BOOST_LOG_TRIVIAL(debug) << "Time elapsed (CGAL): " << elapsed << " ms";

  return rms;
}

float MeshToMeshDistance(const std::vector<Triangle> &a, const std::vector<Triangle> &b) {
  auto start = meshdist::MilliSecondsSinceEpoch();

  float d1 = meshdist_cgal::MeshToMeshDistanceOneDirection(a, b, SAMPLING_DENSITY);
  float d2 = meshdist_cgal::MeshToMeshDistanceOneDirection(b, a, SAMPLING_DENSITY);

  auto elapsed = meshdist::MilliSecondsSinceEpoch() - start;
  BOOST_LOG_TRIVIAL(debug) << "Time elapsed (MeshToMeshDistance): " << elapsed << " ms";
  BOOST_LOG_TRIVIAL(debug) << d1 << ", " << d2;
  return static_cast<float>((d1 + d2) * 0.5);
}
}
