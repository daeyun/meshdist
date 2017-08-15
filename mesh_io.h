//
// Created by daeyun on 5/14/17.
//

#ifndef MESHDIST_MESH_IO_H
#define MESHDIST_MESH_IO_H

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <boost/log/trivial.hpp>

namespace meshdist {
namespace io {
bool ReadTriangles(const std::string &filename,
                   const std::function<void(const std::array<std::array<float, 3>, 3> &)> &triangle_handler);
}
}

#endif //MESHDIST_MESH_IO_H
