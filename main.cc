#define CGAL_HAS_NO_THREADS 1

#include <iostream>
#include <chrono>
#include <ostream>
#include <iomanip>
#include "meshdist.h"
#include "mesh_io.h"
#include "meshdist_cgal.h"

#include <boost/log/expressions.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>


using namespace meshdist;

int main() {
  Eigen::initParallel();
  boost::log::core::get()->set_filter(
      boost::log::trivial::severity >= boost::log::trivial::debug
  );

  std::map<std::string, std::unique_ptr<std::vector<Triangle>>> gt_meshes;
  std::map<std::string, std::unique_ptr<std::vector<Triangle>>> fssr_meshes;
  const int skip = 300;
  const int repeat = 10;
  std::string experiments[] = {"novelview", "novelmodel", "novelclass"};

  int num_max_threads = omp_get_max_threads();
  BOOST_LOG_TRIVIAL(info) << "max_threads: " << num_max_threads;

  for (int k = 1; k <= 2 * num_max_threads; ++k) {
//  for (int k = num_max_threads; k <= num_max_threads; ++k) {
//  for (int k = 2; k <= num_max_threads; k+=10) {
    omp_set_num_threads(k);
    int num_threads = omp_get_max_threads();
    BOOST_LOG_TRIVIAL(info) << "threads: " << num_threads;

    std::vector<double> cgal_elapsed_list;
    std::vector<double> meshdist_elapsed_list;

    for (int l = 0; l < repeat; ++l) {
      int total_count = 0;
      double cgal_elapsed = 0;
      double meshdist_elapsed = 0;

      for (int j = 0; j < 3; ++j) {
        float sum_meshdist = 0;
        float sum_cgal = 0;
        int count = 0;

        for (int i = 0; i < 600; i += skip) {
          std::string basedir = "/data/daeyun/shrec12_recon_links/";

          std::ostringstream dir_stream;
          dir_stream << basedir;
          dir_stream << experiments[j] << "/";
          dir_stream << std::setw(3) << std::setfill('0') << i;

          std::string gt_filename = dir_stream.str() + "/gt.off";
          std::string fssr_filename = dir_stream.str() + "/fssr.off";

          if (gt_meshes.find(gt_filename) == gt_meshes.end()) {
            auto mesh = std::make_unique<std::vector<Triangle>>();
            io::ReadTriangles(gt_filename, [&](auto v) {
              mesh->emplace_back(
                  Vec3(v[0][0], v[0][1], v[0][2]),
                  Vec3(v[1][0], v[1][1], v[1][2]),
                  Vec3(v[2][0], v[2][1], v[2][2])
              );
            });
            gt_meshes[gt_filename] = std::move(mesh);
          }
          std::vector<Triangle> *triangles = gt_meshes[gt_filename].get();

          if (fssr_meshes.find(fssr_filename) == fssr_meshes.end()) {
            auto mesh = std::make_unique<std::vector<Triangle>>();
            io::ReadTriangles(fssr_filename, [&](auto v) {
              mesh->emplace_back(
                  Vec3(v[0][0], v[0][1], v[0][2]),
                  Vec3(v[1][0], v[1][1], v[1][2]),
                  Vec3(v[2][0], v[2][1], v[2][2])
              );
            });
            fssr_meshes[fssr_filename] = std::move(mesh);
          }
          std::vector<Triangle> *triangles2 = fssr_meshes[fssr_filename].get();

          Points3 points, points2;
          meshdist::SamplePointsOnTriangles(*triangles, 900, &points);
          meshdist::SamplePointsOnTriangles(*triangles, 900, &points2);

          float mean_point_distance = (points - points2).colwise().norm().mean();
          BOOST_LOG_TRIVIAL(debug) << "Mean point distance: " << mean_point_distance;

          long start, end;

          start = meshdist::MilliSecondsSinceEpoch();
          float dist_cgal = meshdist_cgal::MeshToMeshDistance(*triangles, *triangles2);
          end = meshdist::MilliSecondsSinceEpoch();
          sum_cgal += dist_cgal / mean_point_distance;
          BOOST_LOG_TRIVIAL(debug) << "cgal: " << dist_cgal;
          cgal_elapsed += end - start;

          start = meshdist::MilliSecondsSinceEpoch();
          // TODO
//          float dist_meshdist = meshdist::MeshToMeshDistance(*triangles, *triangles2);
          float dist_meshdist = 0;
          end = meshdist::MilliSecondsSinceEpoch();
          sum_meshdist += dist_meshdist / mean_point_distance;
          BOOST_LOG_TRIVIAL(debug) << "meshdist: " << dist_meshdist;
          meshdist_elapsed += end - start;

          count++;
          total_count++;
        }

        BOOST_LOG_TRIVIAL(debug) << "MeshToMeshDistance(cgal): " << sum_cgal / count;
        BOOST_LOG_TRIVIAL(debug) << "MeshToMeshDistance(meshdist): " << sum_meshdist / count;
      }

      cgal_elapsed_list.push_back(cgal_elapsed / total_count);
      meshdist_elapsed_list.push_back(meshdist_elapsed / total_count);
    }

    double cgal_best = *std::min_element(cgal_elapsed_list.begin(), cgal_elapsed_list.end());
    double meshdist_best = *std::min_element(meshdist_elapsed_list.begin(), meshdist_elapsed_list.end());

    BOOST_LOG_TRIVIAL(info) << "Elapsed(cgal): " << cgal_best;
    BOOST_LOG_TRIVIAL(info) << "Elapsed(meshdist): " << meshdist_best;

  }

  return 0;
}