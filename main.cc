#include <iostream>
#include <chrono>
#include "meshdist.h"
#include "TriMesh.h"

using namespace meshdist;
using namespace trimesh;

void ReadMesh(const std::string &filename, std::vector<Triangle> *triangles) {
  TriMesh *m = TriMesh::read(filename);
  assert(m != nullptr);
  m->need_faces();
  triangles->reserve(triangles->size() + m->faces.size());
  LOG(m->faces.size());
  for (int j = 0; j < m->faces.size(); ++j) {
    triangles->emplace_back(
        Vec3{m->vertices[m->faces[j][0]][0], m->vertices[m->faces[j][0]][1],
             m->vertices[m->faces[j][0]][2]},
        Vec3{m->vertices[m->faces[j][1]][0], m->vertices[m->faces[j][1]][1],
             m->vertices[m->faces[j][1]][2]},
        Vec3{m->vertices[m->faces[j][2]][0], m->vertices[m->faces[j][2]][1],
             m->vertices[m->faces[j][2]][2]}
    );
  }
  delete m;
}

int main() {
  Eigen::initParallel();

//  std::string filename =
//      "/media/daeyun/Data/data/shapenetcore/ShapeNetCore.v2/03624134/3b1f7f066991f2d45969e7cd6a0b6a55/models/model_normalized.obj";
//  std::string filename2 =
//      "/media/daeyun/Data/data/shapenetcore/ShapeNetCore.v2/03624134/3cbec0d04a115d9c239ffdb3f2fb535d/models/model_normalized.obj";
  std::string filename = "/media/daeyun/Data/data/shapenetcore/ShapeNetCore.v2/03001627/1a6f615e8b1b5ae4dbbc9440457e303e/models/model_normalized.obj";
  std::string filename2 = "/media/daeyun/Data/data/shapenetcore/ShapeNetCore.v2/03001627/1a8bbf2994788e2743e99e0cae970928/models/model_normalized.obj";

  std::vector<Triangle> triangles;
  ReadMesh(filename, &triangles);

  std::vector<Triangle> triangles2;
  ReadMesh(filename2, &triangles2);

  auto start = MilliSecondsSinceEpoch();

  float d1 = MeshToMeshDistance(triangles, triangles2);

  std::cout << d1 << std::endl;

  auto end = MilliSecondsSinceEpoch() - start;
  LOG(end);

  return 0;

#if 0
  std::vector<Triangle> triangles;
  int num_triangles = 2000;
  int num_points = 2000;
  for (int i = 0; i < num_triangles; ++i) {
    triangles.emplace_back(meshdist::Vec3::Random(),
                           Vec3::Random(),
                           Vec3::Random());
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
#endif
}