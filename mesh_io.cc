//
// Created by daeyun on 5/14/17.
//
#include "mesh_io.h"


namespace meshdist {
namespace io {
bool ReadTriangles(const std::string &filename,
                   const std::function<void(const std::array<std::array<float, 3>, 3> &)> &triangle_handler) {
  Assimp::Importer importer;
  const aiScene *scene = importer.ReadFile(filename, aiProcess_Triangulate);
  BOOST_LOG_TRIVIAL(debug) << "Importing " << filename;

  if (!scene) {
    BOOST_LOG_TRIVIAL(error) << importer.GetErrorString();
    return false;
  }

  for (int i = 0; i < scene->mNumMeshes; ++i) {
    const aiMesh *mesh = scene->mMeshes[i];
    for (int j = 0; j < mesh->mNumFaces; ++j) {
      auto face = mesh->mFaces[j];
      if (face.mNumIndices == 3) {
        for (int k = 0; k < 3; ++k) {
          if (face.mIndices[k] >= mesh->mNumVertices) {
            BOOST_LOG_TRIVIAL(warning) << "Invalid vertex index found. Skipping.";
            continue;
          }
        }
        auto a = mesh->mVertices[face.mIndices[0]];
        auto b = mesh->mVertices[face.mIndices[1]];
        auto c = mesh->mVertices[face.mIndices[2]];
        triangle_handler({std::array<float, 3>{a.x, a.y, a.z},
                          std::array<float, 3>{b.x, b.y, b.z},
                          std::array<float, 3>{c.x, c.y, c.z}});
      }
    }
  }

  return true;
}
}
}
