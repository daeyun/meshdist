#include <iostream>
#include <chrono>
#include "meshdist.h"

long MilliSecondsSinceEpoch() {
  return std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

int main() {

  std::vector<Triangle> triangles;
  int num_triangles = 2000;
  int num_points = 2000;
  for (int i = 0; i < num_triangles; ++i) {
    triangles.emplace_back(Vec3::Random(), Vec3::Random(), Vec3::Random());
  }
  auto points = Points3d::Random(3, num_points);

  std::vector<std::vector<Float>> all_distances;

  long best = std::numeric_limits<long>::max();
  for (int i = 0; i < 100; ++i) {
    auto start = MilliSecondsSinceEpoch();
    DistanceToTriangles(triangles, points, &all_distances);
    auto end = MilliSecondsSinceEpoch() - start;
    if (end < best) {
      best = end;
    }
  }
  LOG(best);

  return 0;
}