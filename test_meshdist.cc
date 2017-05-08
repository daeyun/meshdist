//
// Created by daeyun on 4/28/17.
//

#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"

#include "meshdist.h"

TEST_CASE("Rotate around x axis") {
  auto R = AxisRotation(M_PI, {1, 0, 0});
  Vec3 p{0, 1, 1};
  auto Rp = R * p;
  REQUIRE(Rp.isApprox(Vec3(0, -1, -1)));
}

TEST_CASE("Rotation matrix identity") {
  auto R = AxisRotation(M_PI, {1, -2, 3});
  Vec3 p = Vec3::Random();
  auto Rp = R * (R * p);
  REQUIRE(Rp.isApprox(p));
}

TEST_CASE("Check for zero vector") {
  REQUIRE_THROWS(AxisRotation(M_PI, {0, 0, 0}));
  REQUIRE_THROWS(AxisRotation(M_PI, {0, 0, 1e-12}));
}

TEST_CASE("yz plane alignment identity") {
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  Triangle triangle_copy = triangle;
  Mat34 M;
  YZTransform(&triangle, &M);
  REQUIRE(triangle_copy.a.isApprox(triangle.a));
  REQUIRE(triangle_copy.b.isApprox(triangle.b));
  REQUIRE(triangle_copy.c.isApprox(triangle.c));
}

TEST_CASE("yz plane alignment changes value") {
  Triangle triangle({0, 0, 1}, {0, 0, 2}, {0, 3, 0});
  Triangle triangle_copy = triangle;
  Mat34 M;
  YZTransform(&triangle, &M);
  REQUIRE(!triangle_copy.a.isApprox(triangle.a));
  REQUIRE(!triangle_copy.b.isApprox(triangle.b));
  REQUIRE(!triangle_copy.c.isApprox(triangle.c));
}

TEST_CASE("yz plane alignment inverse") {
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  Triangle triangle_copy = triangle;

  for (int i = 0; i < 10; ++i) {
    Vec3 axis;
    axis.setRandom();

    Mat33 R = AxisRotation(
        (std::rand() / static_cast<Float>(RAND_MAX)) * M_PI * 2, axis);
    triangle.ApplyR(R);
    triangle.Translate(Vec3::Random());
    Triangle transformed = triangle;

    Mat34 M;
    YZTransform(&triangle, &M);
    REQUIRE(triangle_copy.a.isApprox(triangle.a));
    REQUIRE(triangle_copy.b.isApprox(triangle.b));
    REQUIRE(triangle_copy.c.isApprox(triangle.c));

    transformed.ApplyRt(M);
    REQUIRE(transformed.a.isApprox(triangle.a));
    REQUIRE(transformed.b.isApprox(triangle.b));
    REQUIRE(transformed.c.isApprox(triangle.c));
  }
}

TEST_CASE("closest triangle part is A") {
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  Points3d points = -Points3d::Random(3, 200).cwiseAbs();
  auto parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::A);
  }
}

TEST_CASE("closest triangle part is AB") {
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  Points3d points = Points3d::Random(3, 200).cwiseAbs();
  points.row(1) *= -1;
  points.row(2) *= 2;
  auto parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::AB);
  }
}

TEST_CASE("closest triangle part is CA") {
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  Points3d points = Points3d::Random(3, 200).cwiseAbs();
  points.row(1) *= 3;
  points.row(2) *= -1;
  auto parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::CA);
  }
}

TEST_CASE("closest triangle part is BC") {
  Triangle triangle({0, 0, 0}, {0, 0, 1}, {0, 1, 0});
  Points3d points = Points3d::Random(3, 200).cwiseAbs();
  for (int i = 0; i < points.cols(); ++i) {
    if (points.col(i).bottomRows(2).sum() < 1) {
      points.col(i).bottomRows(2) = 1 - points.col(i).bottomRows(2).array();
    }
  }
  auto parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::BC);
  }
}

TEST_CASE("closest triangle part is ABC") {
  Triangle triangle({0, 0, 0}, {0, 0, 1}, {0, 1, 0});
  Points3d points = Points3d::Random(3, 200).cwiseAbs();
  for (int i = 0; i < points.cols(); ++i) {
    if (points.col(i).bottomRows(2).sum() > 1) {
      points.col(i).bottomRows(2) = 1 - points.col(i).bottomRows(2).array();
    }
  }
  auto parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::ABC);
  }
}

TEST_CASE("closest triangle part is B") {
  Triangle triangle({0, 0, 0}, {0, 0, 1}, {0, 1, 0});
  Points3d points = Points3d::Random(3, 200).cwiseAbs();
  for (int i = 0; i < points.cols(); ++i) {
    if (points.col(i).bottomRows(2).sum() > 1) {
      points.col(i).bottomRows(2) = 1 - points.col(i).bottomRows(2).array();
    }
  }
  points.row(2) = 2 - points.row(2).array();
  auto parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::B);
  }
  points = Points3d::Random(3, 200).cwiseAbs();
  points.row(2).array() += 1;
  points.row(1).array() *= -1;
  parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::B);
  }
}

TEST_CASE("closest triangle part is C") {
  Triangle triangle({0, 0, 0}, {0, 0, 1}, {0, 1, 0});
  Points3d points = Points3d::Random(3, 200).cwiseAbs();
  for (int i = 0; i < points.cols(); ++i) {
    if (points.col(i).bottomRows(2).sum() > 1) {
      points.col(i).bottomRows(2) = 1 - points.col(i).bottomRows(2).array();
    }
  }
  points.row(1) = 2 - points.row(1).array();
  auto parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::C);
  }
  points = Points3d::Random(3, 200).cwiseAbs();
  points.row(1).array() += 1;
  points.row(2).array() *= -1;
  parts = FindClosestPart(triangle, points);
  for (auto &&part : parts) {
    REQUIRE(part == Part::C);
  }
}

TEST_CASE("closest point on triangle") {
  std::vector<Triangle> triangles;
  for (int i = 0; i < 10; ++i) {
    triangles.emplace_back(Vec3::Random(), Vec3::Random(), Vec3::Random());
  }
  Points3d points = Points3d::Random(3, 10);

  std::vector<std::vector<Float>> distances;
  DistanceToTriangles(triangles, points, &distances);

  std::vector<std::vector<Float>> sample_distances;
  std::vector<std::vector<Part>> parts;
  for (auto &&triangle : triangles) {
    std::vector<Vec3> triangle_points;
    for (int i = 0; i < 10000; ++i) {
      triangle_points.push_back(triangle.SamplePoint());
    }
    std::vector<Float> sd;
    for (int j = 0; j < points.cols(); ++j) {
      Float closest = std::numeric_limits<Float>::max();
      for (auto &&triangle_point : triangle_points) {
        Float d = (triangle_point - points.col(j)).norm();
        if (d < closest) {
          closest = d;
        }
      }
      sd.push_back(closest);
    }
    sample_distances.push_back(sd);

    Mat34 T;
    Triangle verts = triangle;
    YZTransform(&verts, &T);

    Points3d Tpoints = ApplyMat34(T, points);
    parts.push_back(FindClosestPart(verts, Tpoints));
  }

  constexpr Float kInf = std::numeric_limits<Float>::max();
  std::map<Part, Float> max_diff{
      {Part::A, -kInf}, {Part::B, -kInf}, {Part::C, -kInf},
      {Part::AB, -kInf}, {Part::BC, -kInf}, {Part::CA, -kInf},
      {Part::ABC, -kInf},
  };
  std::map<Part, int> counts;

  for (int k = 0; k < distances.size(); ++k) {
    for (int i = 0; i < distances[k].size(); ++i) {
      Float diff = std::abs(sample_distances[k][i] - distances[k][i]);
      if (diff > max_diff[parts[k][i]]) {
        max_diff[parts[k][i]] = diff;
      }
      counts[parts[k][i]]++;
    }
  }

  // Make sure there's more than one example of each type.
  REQUIRE(counts[Part::A] > 1);
  REQUIRE(counts[Part::B] > 1);
  REQUIRE(counts[Part::C] > 1);
  REQUIRE(counts[Part::AB] > 1);
  REQUIRE(counts[Part::BC] > 1);
  REQUIRE(counts[Part::CA] > 1);
  REQUIRE(counts[Part::ABC] > 1);

  REQUIRE(max_diff[Part::A] < 0.03);
  REQUIRE(max_diff[Part::B] < 0.03);
  REQUIRE(max_diff[Part::C] < 0.03);
  REQUIRE(max_diff[Part::AB] < 0.03);
  REQUIRE(max_diff[Part::BC] < 0.03);
  REQUIRE(max_diff[Part::CA] < 0.03);
  REQUIRE(max_diff[Part::ABC] < 0.03);
}

TEST_CASE("closest point to BC") {
  Triangle triangle({0, 0, 0}, {0, 0, 1}, {0, 1, 0});
  Vec3 point{0, 1 - 0.1, 1 + 0.1};

  auto expected = Approx(std::sqrt(2.0) / 2);

  std::vector<std::vector<Float>> distances;
  DistanceToTriangles({triangle}, Points3d{point}, &distances);
  REQUIRE(distances[0][0] == expected);

  for (int i = 0; i < 20; ++i) {
    Mat33 R = AxisRotation(
        (std::rand() / static_cast<Float>(RAND_MAX)) * M_PI * 2,
        Vec3::Random());
    Vec3 t = Vec3::Random();

    triangle.ApplyR(R);
    triangle.Translate(t);
    point = (R * point + t).eval();

    distances.clear();
    DistanceToTriangles({triangle}, Points3d{point}, &distances);
    REQUIRE(distances[0][0] == expected);
  }
}

TEST_CASE("closest point to CA") {
  Triangle triangle({0, 0, 0}, {0, 0, 1}, {0, 1, 1});
  Vec3 point{0, 1 + 0.1, 0.1};

  auto expected = Approx(std::sqrt(2.0) / 2);

  std::vector<std::vector<Float>> distances;
  DistanceToTriangles({triangle}, Points3d{point}, &distances);
  REQUIRE(distances[0][0] == expected);

  for (int i = 0; i < 20; ++i) {
    Mat33 R = AxisRotation(
        (std::rand() / static_cast<Float>(RAND_MAX)) * M_PI * 2,
        Vec3::Random());
    Vec3 t = Vec3::Random();

    triangle.ApplyR(R);
    triangle.Translate(t);
    point = (R * point + t).eval();

    distances.clear();
    DistanceToTriangles({triangle}, Points3d{point}, &distances);
    REQUIRE(distances[0][0] == expected);
  }
}
