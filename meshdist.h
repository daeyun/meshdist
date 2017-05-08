//
// Created by daeyun on 4/28/17.
//

#ifndef MESHDIST_MESHDIST_H
#define MESHDIST_MESHDIST_H

#include <random>
#include <Eigen/Dense>
#include <boost/optional.hpp>

typedef float Float;

using Eigen::Dynamic;
using Vec = Eigen::Matrix<Float, Dynamic, 1>;
using Vec2 = Eigen::Matrix<Float, 2, 1>;
using Vec3 = Eigen::Matrix<Float, 3, 1>;
using Vec4 = Eigen::Matrix<Float, 4, 1>;
using Mat44 = Eigen::Matrix<Float, 4, 4>;
using Mat34 = Eigen::Matrix<Float, 3, 4>;
using Mat33 = Eigen::Matrix<Float, 3, 3>;
using Points3d = Eigen::Matrix<Float, 3, Dynamic>;
using Points2d = Eigen::Matrix<Float, 2, Dynamic>;
using Points2i = Eigen::Matrix<int, 2, Dynamic>;

#define LOG(msg) std::cout << (#msg) << ": " << (msg) << std::endl;

constexpr Float kEpsilon = 1e-9;

enum class Part { A, B, C, AB, BC, CA, ABC };

auto &RandomEngine() {
  thread_local static std::mt19937 engine{std::random_device{}()};
  return engine;
}

auto Random() {
  thread_local static std::uniform_real_distribution<Float> dist{0, 1};
  return dist(RandomEngine());
}

class Triangle {
 public:
  Triangle(const Vec3 &a, const Vec3 &b, const Vec3 &c) : a(a), b(b), c(c) {}

  void ApplyRt(const Mat34 &M) {
    auto R = M.leftCols<3>();
    auto t = M.col(3);
    a = (R * a + t).eval();
    b = (R * b + t).eval();
    c = (R * c + t).eval();
    ClearCache();
  }

  void ApplyR(const Mat33 &M) {
    a = (M * a).eval();
    b = (M * b).eval();
    c = (M * c).eval();
    ClearCache();
  }

  void Translate(const Vec3 &dxyz) {
    a += dxyz;
    b += dxyz;
    c += dxyz;
    ClearCache();
  }

  void Print() const {
    std::cout << a.transpose() << ",      "
              << b.transpose() << ",      "
              << c.transpose() << std::endl;
  }

  Vec3 SamplePoint() const {
    Vec2 v12{Random(), Random()};
    if (v12.sum() > 1) {
      v12 = 1 - v12.array();
    }
    return a + (v12(0) * ab()) + (v12(1) * ac());
  }

  Float Area() const {
    return ab().cross(ac()).norm() * 0.5f;
  }

  void ClearCache() {
    ab_.reset();
    ac_.reset();
    sin_cos_bc_.reset();
    sin_cos_ca_.reset();
  }

  const Vec3 &ab() const {
    if (!ab_) { ab_ = b - a; }
    return *ab_;
  }

  const Vec3 &ac() const {
    if (!ac_) { ac_ = c - a; }
    return *ac_;
  }

  const std::pair<Float, Float> &sin_cos_bc() const {
    if (!sin_cos_bc_) {
      Float theta_bc = -std::atan2(c(2) - b(2), c(1));
      sin_cos_bc_.emplace(std::sin(theta_bc), std::cos(theta_bc));
    }
    return *sin_cos_bc_;
  };

  const std::pair<Float, Float> &sin_cos_ca() const {
    if (!sin_cos_ca_) {
      Float theta_ca = -std::atan2(c(2), c(1));
      sin_cos_ca_.emplace(std::sin(theta_ca), std::cos(theta_ca));
    }
    return *sin_cos_ca_;
  };

  Vec3 a, b, c;

 private:
  mutable boost::optional<Vec3> ab_, ac_;
  mutable boost::optional<std::pair<Float, Float>> sin_cos_bc_, sin_cos_ca_;
};

Points3d ApplyMat34(const Mat34 &M, const Points3d &points) {
  return (M.leftCols<3>() * points).colwise() + M.col(3);
}

Mat33 AxisRotation(Float angle, const Vec3 &axis) {
  if (axis.isZero()) {
    throw std::invalid_argument("Zero vector");
  }
  if (std::abs(angle) < kEpsilon) {
    return Mat33::Identity();
  }
  Vec3 rotation_axis = axis.normalized();
  Mat33 R = Mat33::Zero();
  Float cosa = std::cos(angle);
  Float sina = std::sin(angle);
  R.diagonal().fill(cosa);
  R += (rotation_axis * rotation_axis.transpose()) * (1.0 - cosa);
  Vec3 v = rotation_axis * sina;
  R(0, 1) += -v[2];
  R(0, 2) += v[1];
  R(1, 0) += v[2];
  R(1, 2) += -v[0];
  R(2, 0) += -v[1];
  R(2, 1) += v[0];
  return R;
}

Mat33 ZRotation(Float angle) {
  if (std::abs(angle) < kEpsilon) {
    return Mat33::Identity();
  }
  Mat33 R = Mat33::Zero();
  Float cosa = std::cos(angle);
  Float sina = std::sin(angle);
  R(0, 0) = cosa;
  R(0, 1) = -sina;
  R(1, 0) = sina;
  R(1, 1) = cosa;
  R(2, 2) = 1;
  return R;
}

void YZTransform(Triangle *verts, Mat34 *M = nullptr) {
  // 1. Translate so that A = (0, 0, 0).
  Vec3 translation = verts->a;
  verts->b.array() -= translation.array();
  verts->c.array() -= translation.array();
  verts->a = {0, 0, 0};

  // 2. Rotate so that AB aligns with the positive z axis.
  auto z = Vec3(0, 0, 1);

  // Angle between AB and z-axis.
  Vec3 ab = verts->b.normalized();
  Vec3 rot_axis = ab.cross(z);
  // More numerically stable than arccos.
  Float cos_theta = ab.dot(z);
  Float theta = std::atan2(rot_axis.norm(), cos_theta);

  // If ab is already on the z-axis, use y-axis.
  if (ab.head<2>().isZero()) {
    rot_axis = {0, 1, 0};
  }

  Mat33 R1 = AxisRotation(theta, rot_axis);
  Vec2 cxy = R1.topLeftCorner<2, 3>() * verts->c;

  // 3. Rotate around z-axis so that the triangle is on the yz plane.
  Float theta_z = static_cast<Float>(M_PI_2) - std::atan2(cxy(1), cxy(0));

  Mat33 R2 = ZRotation(theta_z);
  Mat33 R = R2 * R1;

  verts->ApplyR(R);

  if (M) {
    M->topLeftCorner<3, 3>() = R;
    M->col(3).array() = R * -translation;
  }
}

bool IsLeft(const Vec2 &p, const Vec2 &origin, const Vec2 &direction) {
  auto po = p - origin;  //  `auto` is important for lazy evaluation.
  auto d0 = po(0) * direction(1);
  auto d1 = po(1) * direction(0);
  return d0 - d1 < 0.0;
}

Vec2 CCRotate2d(const Vec2 &v) {
  return {-v[1], v[0]};
}

std::vector<Part> FindClosestPart(Triangle &yz_triangle,
                                  const Points3d &points) {
  auto p1 = yz_triangle.a.tail<2>();
  auto p2 = yz_triangle.b.tail<2>();
  auto p3 = yz_triangle.c.tail<2>();

  Vec2 d31 = p1 - p3;
  Vec2 d31_cc = CCRotate2d(d31);
  Vec2 d12 = p2 - p1;
  Vec2 d12_cc = CCRotate2d(d12);
  Vec2 d23 = p3 - p2;
  Vec2 d23_cc = CCRotate2d(d23);

  std::vector<Part> parts;
  const size_t n = static_cast<size_t>(points.cols());
  parts.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    Vec2 point_yz = points.col(i).tail<2>();

    bool s2 = IsLeft(point_yz, p1, d31_cc);
    bool s4 = IsLeft(point_yz, p1, d12_cc);
    if (!s2 && s4) {
      parts.push_back(Part::A);
      continue;
    }

    bool s5 = IsLeft(point_yz, p2, d12_cc);
    bool s7 = IsLeft(point_yz, p2, d23_cc);
    if (!s5 && s7) {
      parts.push_back(Part::B);
      continue;
    }

    bool s1 = IsLeft(point_yz, p3, d31_cc);
    bool s8 = IsLeft(point_yz, p3, d23_cc);
    if (!s8 && s1) {
      parts.push_back(Part::C);
      continue;
    }

    bool s3 = IsLeft(point_yz, p1, d12);
    if (s3 && !s4 && s5) {
      parts.push_back(Part::AB);
      continue;
    }

    bool s6 = IsLeft(point_yz, p2, d23);
    if (s6 && !s7 && s8) {
      parts.push_back(Part::BC);
      continue;
    }

    bool s0 = IsLeft(point_yz, p3, d31);
    if (s0 && !s1 && s2) {
      parts.push_back(Part::CA);
      continue;
    }

    parts.push_back(Part::ABC);
  }

  return parts;
}

void DistanceToTriangles(const std::vector<Triangle> &triangles,
                         const Points3d &points,
                         std::vector<std::vector<Float>> *all_distances) {
  for (auto &&triangle : triangles) {
    Mat34 T;
    Triangle verts = triangle;
    YZTransform(&verts, &T);

    Points3d Tpoints = ApplyMat34(T, points);
    auto parts = FindClosestPart(verts, Tpoints);

    std::vector<Float> distances;
    distances.reserve(parts.size());
    for (int i = 0; i < parts.size(); ++i) {
      Part part = parts[i];
      Vec3 Tpoint = Tpoints.col(i).eval();

      Float dist;
      switch (part) {
        case Part::A:
          dist = Tpoint.norm();
          break;
        case Part::B:
          Tpoint(2) -= verts.b(2);
          dist = Tpoint.norm();
          break;
        case Part::C:
          Tpoint.bottomRows(2) -= verts.c.bottomRows(2);
          dist = Tpoint.norm();
          break;
        case Part::AB:
          dist = Tpoint.topRows(2).norm();
          break;
        case Part::BC: {
          const auto sin_cos = verts.sin_cos_bc();
          Float z;
          z = sin_cos.first * Tpoint(1)
              + sin_cos.second * (Tpoint(2) - verts.b(2));
          dist = std::sqrt(z * z + Tpoint(0) * Tpoint(0));
          break;
        }
        case Part::CA: {
          const auto sin_cos = verts.sin_cos_ca();
          Float z;
          z = sin_cos.first * Tpoint(1) + sin_cos.second * Tpoint(2);
          dist = std::sqrt(z * z + Tpoint(0) * Tpoint(0));
          break;
        }
        case Part::ABC:
          dist = std::abs(Tpoint(0));
          break;
      }
      distances.push_back(dist);
    }
    all_distances->push_back(distances);
  }
}

void SamplePointsOnTriangles(const std::vector<Triangle> &triangles,
                             int num_samples, Points3d *points) {
  std::vector<Float> areas;
  areas.reserve(triangles.size());
  for (auto &&triangle : triangles) {
    areas.push_back(triangle.Area());
  }
  std::vector<int> counts(triangles.size());
  std::discrete_distribution<> distribution(std::begin(areas), std::end(areas));

  for (int i = 0; i < num_samples; ++i) {
    ++counts[distribution(RandomEngine())];
  }

  points->resize(3, num_samples);
  int k = 0;
  for (int i = 0; i < triangles.size(); ++i) {
    for (int j = 0; j < counts[i]; ++j) {
      points->col(k) = triangles[i].SamplePoint();
      ++k;
    }
  }
}

#endif //MESHDIST_MESHDIST_H
