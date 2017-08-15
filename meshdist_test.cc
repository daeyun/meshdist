//
// Created by daeyun on 4/28/17.
//

#include "gtest/gtest.h"

#include "meshdist.h"

using namespace meshdist;

double RandomAngle() {
  return (std::rand() / static_cast<double>(RAND_MAX)) * 2 * M_PI - M_PI;
}

Mat34 RandomRt() {
  return AxisRotation(RandomAngle(), Vec3::Random(), Vec3::Random());
}

TEST(YZTransform, Identity) {
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  Triangle triangle_copy = triangle;
  Mat34 M;
  YZTransform(&triangle, &M);
  for (int i = 0; i < 3; ++i) {
    ASSERT_TRUE(triangle_copy[i].isApprox(triangle[i]));
  }
}

TEST(YZTransform, Alignment) {
  Triangle triangle({0, 0, 1}, {0, 0, 2}, {0, 3, 0});
  Triangle triangle_copy = triangle;
  Mat34 M;
  YZTransform(&triangle, &M);
  for (int i = 0; i < 3; ++i) {
    auto diff = (triangle_copy[i] - triangle[i]).cwiseAbs().maxCoeff();
    ASSERT_GT(diff, 0.1);
  }
}

TEST(YZTransform, FlipZ) {
  Triangle triangle({0, 0, 0}, {0, 0, -2}, {0, 3, 0});
  Triangle triangle_copy = triangle;
  Mat34 M;
  YZTransform(&triangle, &M);
  ASSERT_TRUE(triangle.b.isApprox(Vec3(0, 0, 2)));
}

TEST(YZTransform, PositiveAxes) {
  for (int i = 0; i < 500; ++i) {
    Triangle triangle(Vec3::Random(), Vec3::Random(), Vec3::Random());
    Mat34 M;
    YZTransform(&triangle, &M);
    ASSERT_FLOAT_EQ(triangle.a[0], 0);
    ASSERT_FLOAT_EQ(triangle.a[1], 0);
    ASSERT_FLOAT_EQ(triangle.a[2], 0);

    ASSERT_NEAR(triangle.b[0], 0, 1e-6);
    ASSERT_NEAR(triangle.b[1], 0, 1e-6);
    ASSERT_GT(triangle.b[2], 1e-6);

    ASSERT_NEAR(triangle.c[0], 0, 1e-6);
    ASSERT_GT(triangle.c[1], 1e-6);
    // vertex c can be in the negative z axis.
  }
}

TEST(YZTransform, Inverse) {
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  Triangle triangle_copy = triangle;

  for (int i = 0; i < 500; ++i) {
    Vec3 axis;
    axis.setRandom();

    Mat34 rand_Rt = RandomRt();
    triangle.ApplyRt(rand_Rt);
    Triangle transformed = triangle;

    Mat34 M;
    YZTransform(&triangle, &M);

    transformed.ApplyRt(M);

    for (int i = 0; i < 3; ++i) {
      double diff;

      diff = (triangle_copy[i] - triangle[i]).cwiseAbs().maxCoeff();
      ASSERT_LT(diff, 1e-4);

      diff = (transformed[i] - triangle[i]).cwiseAbs().maxCoeff();
      ASSERT_LT(diff, 1e-6);
    }

  }
}

TEST(YZTransform, Illconditioned) {
  Triangle triangle1({0, 0, 0}, {0, 0, 0}, {0, 0, 0});
  Triangle triangle2({0, 0, 0}, {0, 0, 0}, {1, 0, 0});
  Triangle triangle3({0, 0, 0}, {0, 0, 0}, {0, 1, 0});
  Triangle triangle4({0, 0, 0}, {0, 0, 0}, {0, 0, 1});
  Triangle triangle5({0, 0, 0}, {0, 0, 0}, {-1, 0, 0});
  Triangle triangle6({0, 0, 0}, {0, 0, 0}, {0, -1, 0});
  Triangle triangle7({0, 0, 0}, {0, 0, 0}, {0, 0, -1});
  Triangle triangle8({0, 0, 0}, {2, 0, 0}, {-1, 0, 0});
  Triangle triangle9({0, 0, 0}, {0, 2, 2}, {0, 3, 3});

  EXPECT_TRUE(triangle1.SamplePoint().isZero());

  Mat34 M;

  YZTransform(&triangle1, &M);
  EXPECT_TRUE(M.allFinite());

  YZTransform(&triangle2, &M);
  EXPECT_TRUE(M.allFinite());

  YZTransform(&triangle3, &M);
  EXPECT_TRUE(M.allFinite());

  YZTransform(&triangle4, &M);
  EXPECT_TRUE(M.allFinite());

  YZTransform(&triangle5, &M);
  EXPECT_TRUE(M.allFinite());

  YZTransform(&triangle6, &M);
  EXPECT_TRUE(M.allFinite());

  YZTransform(&triangle7, &M);
  EXPECT_TRUE(M.allFinite());

  YZTransform(&triangle8, &M);
  EXPECT_TRUE(M.allFinite());

  YZTransform(&triangle9, &M);
  EXPECT_TRUE(M.allFinite());
}

TEST(IntegrationTest, Sampling) {
  Triangle triangle(Vec3::Random(), Vec3::Random(), Vec3::Random());
  Points3 points;
  SamplePointsOnTriangles({triangle}, 1000, &points);

  Vec dist;
  MinimumDistanceToTriangles({triangle}, points, &dist);

}

TEST(DistanceToLineTest, Simple) {
  Vec3 A{0, 0, 1};
  Vec3 B{0, 1, 0};

  Vec3 u = (B - A).normalized();
  Vec3 P{0, 0.0001, 0.0001};

  Mat34 Rt = AxisRotation(M_PI_2, A, u);
  Float expected = DistanceToLine(A, u, P);

  P = ApplyMat34(Rt, P).eval();

  Float out = DistanceToLine(A, u, P);
  EXPECT_FLOAT_EQ(expected, out);
}

TEST(DistanceToLineTest, RotateAroundZAxis) {
  for (int i = 0; i < 100; ++i) {
    Vec3 A{0, 1, 0};
    Vec3 u{0, 0, 1};
    Vec3 P{0, 0.1, 0};

    Float expected = (P - A).head<2>().norm();

    // Rotate and check the distance to the axis stays the same.
    Mat34 Rt = AxisRotation(RandomAngle(), A, u);
    P = ApplyMat34(Rt, P).eval();

    Float out = DistanceToLine(A, u, P);
    ASSERT_FLOAT_EQ(out, expected);
  }
}

TEST(DistanceToLineTest, RotateAroundRandomAxis) {
  for (int i = 0; i < 100; ++i) {
    Vec3 A = Vec3::Random();
    Vec3 u = Vec3::Random();
    Vec3 P = Vec3::Random();

    Float expected = DistanceToLine(A, u, P);

    // Rotate and check the distance to the axis stays the same.
    Mat34 Rt = AxisRotation(RandomAngle(), A, u);
    P = ApplyMat34(Rt, P).eval();

    Float out = DistanceToLine(A, u, P);
    ASSERT_NEAR(out, expected, 1e-6);
  }
}

TEST(RotationAround3dLine, RotateBy180) {
  Vec3 b{0, 0, 1};
  Vec3 c{0, std::sqrt(3), 0};
  Vec3 point{0, 0, 0};

  Mat34 Rt = AxisRotation(M_PI, b, c - b);
  Vec3 transformed = ApplyMat34(Rt, point);

  Mat22 M = Eigen::Rotation2D<Float>(DegreesToRadians(120.0)).matrix();
  Vec2 yz = M * (point.tail<2>() - b.tail<2>()) + b.tail<2>();

  EXPECT_NEAR(transformed(0), 0, 1e-7);
  EXPECT_FLOAT_EQ(transformed(1), yz(0));
  EXPECT_FLOAT_EQ(transformed(2), yz(1));
}

TEST(RotationAround3dLine, RotateBy90) {
  Vec3 b{0, 0, 1};
  Vec3 c{0, std::sqrt(3), 0};
  Vec3 point{0, 0, 0};
  Mat34 Rt = AxisRotation(M_PI_2, b, c - b);
  Vec3 transformed = ApplyMat34(Rt, point);
  Vec2 yz = (b + (c - b).normalized() * 0.5).tail<2>();

  EXPECT_FLOAT_EQ(-std::sqrt(3) * 0.5, transformed(0));
  EXPECT_FLOAT_EQ(yz(0), transformed(1));
  EXPECT_FLOAT_EQ(yz(1), transformed(2));
}

TEST(FindClosestPart, ClosestEdges) {
  // Assume already transformed to fit in the yz plane.
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  // Point in triangle.
  Vec3 point{0, 0.1, 0.1};

  for (int i = 0; i < 1000; ++i) {
    double angle = (i == 0) ? 0 : RandomAngle();

    for (const auto &tup: {std::tuple<Vec3, Vec3, Part>{triangle.a, triangle.b, Part::AB},
                           std::tuple<Vec3, Vec3, Part>{triangle.b, triangle.c, Part::BC},
                           std::tuple<Vec3, Vec3, Part>{triangle.c, triangle.a, Part::CA}}) {
      const Part expected_part = (std::abs(angle) < M_PI_2) ? Part::ABC : std::get<2>(tup);
      Mat34 Rt = AxisRotation(angle, std::get<0>(tup), std::get<1>(tup) - std::get<0>(tup));
      Vec3 rotated_point = ApplyMat34(Rt, Points3{point});

      float pdist = DistanceToLine(std::get<0>(tup), std::get<1>(tup) - std::get<0>(tup), point);
      float out = DistanceToLine(std::get<0>(tup), std::get<1>(tup) - std::get<0>(tup), rotated_point);

      // Sanity check to make sure the transformation is valid.
      EXPECT_NEAR(pdist, out, 1e-6);

      Part part_found = FindClosestPart(triangle, rotated_point)[0];
      ASSERT_EQ(expected_part, part_found);

      Float expected;
      if (expected_part == Part::ABC) {
        expected = std::abs(rotated_point(0));
      } else {
        expected = DistanceToLine(std::get<0>(tup), std::get<1>(tup) - std::get<0>(tup), Points3{rotated_point});
      }

      Mat34 rand_Rt = RandomRt();
      Triangle rand_triangle = triangle;
      rand_triangle.ApplyRt(rand_Rt);
      Vec3 rand_point = ApplyMat34(rand_Rt, rotated_point);

      Vec distances;
      MinimumDistanceToTriangles({rand_triangle}, Points3{rand_point}, &distances);

      ASSERT_EQ(1, distances.cols());
      ASSERT_EQ(1, distances.rows());

      Float out_dist = distances(0);

      ASSERT_NEAR(expected, out_dist, 1e-5);
    }
  }
}

TEST(FindClosestPart, ClosestVertices) {
  // Assume already transformed to fit in the yz plane.
  Triangle triangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  // Point in triangle.
  Vec3 point{0, 0.1, 0.1};
  Vec3 normal{1, 0, 0};

  int count = 0;
  const int n = 2000;

  for (int i = 0; i < n; ++i) {
    double angle = (i == 0) ? 0 : RandomAngle();

    for (const auto &tup: {std::tuple<Vec3, Vec3, Part>{triangle.a, triangle.b, Part::A},
                           std::tuple<Vec3, Vec3, Part>{triangle.b, triangle.c, Part::B},
                           std::tuple<Vec3, Vec3, Part>{triangle.c, triangle.a, Part::C}}) {
      const Part closest_part = (std::abs(angle) < M_PI_2) ? Part::ABC : std::get<2>(tup);
      const auto closest_part_position = std::get<0>(tup);

      Vec3 rotation_axis = (closest_part_position - point).cross(normal);

      Mat34 Rt = AxisRotation(angle, closest_part_position, rotation_axis);
      Vec3 rotated_point = ApplyMat34(Rt, Points3{point});

      float pdist = (closest_part_position - point).norm();
      float out = (closest_part_position - rotated_point).norm();

      // Sanity check to make sure the transformation is valid.
      EXPECT_NEAR(pdist, out, 1e-6);

      Part part_found = FindClosestPart(triangle, rotated_point)[0];
      ASSERT_EQ(closest_part, part_found);

      Float expected_distance;
      if (closest_part == Part::ABC) {
        expected_distance = std::abs(rotated_point(0));
      } else {
        expected_distance = pdist;
        count++;
      }

      Mat34 rand_Rt = RandomRt();
      Triangle rand_triangle = triangle;
      rand_triangle.ApplyRt(rand_Rt);
      Vec3 rand_point = ApplyMat34(rand_Rt, rotated_point);

      Vec distances;
      MinimumDistanceToTriangles({rand_triangle}, Points3{rand_point}, &distances);

      ASSERT_EQ(1, distances.cols());
      ASSERT_EQ(1, distances.rows());

      Float out_dist = distances(0);

      ASSERT_NEAR(expected_distance, out_dist, 1e-5);
    }
  }

  ASSERT_GT(count, n);
}

TEST(SamplePointsOnTriangle, BarycentricArea) {
  for (int j = 0; j < 100; ++j) {
    Triangle triangle(Vec3::Random(), Vec3::Random(), Vec3::Random());
    for (int i = 0; i < 100; ++i) {
      Vec3 point = triangle.SamplePoint();
      Float area_sum = Triangle(triangle.a, triangle.b, point).Area() +
          Triangle(triangle.b, triangle.c, point).Area() +
          Triangle(triangle.c, triangle.a, point).Area();
      ASSERT_NEAR(area_sum, triangle.Area(), 1e-6);
    }
  }
}

TEST(DistanceOneDirection, ParallelTriangles) {
  for (int i = 0; i < 1000; ++i) {
    Triangle triangle1(Vec3::Random(), Vec3::Random(), Vec3::Random());
    Vec3 normal = (triangle1.a - triangle1.b).cross(triangle1.a - triangle1.c).normalized();
    Float d = Random();
    Triangle triangle2 = triangle1;
    triangle2.Translate(normal * d);
    Float rms = MeshToMeshDistanceOneDirection3({triangle1}, {triangle2}, 100);
    ASSERT_NEAR(rms, d, 1e-5);
  }
}

TEST(TriangleObject, DistanceToPoint) {
  for (int i = 0; i < 1000; ++i) {
    Triangle triangle1(Vec3::Random(), Vec3::Random(), Vec3::Random());
    Vec3 p = Vec3::Random();
    float d = triangle1.DistanceTo(p);

    Vec dist;
    MinimumDistanceToTriangles({triangle1}, p, &dist);
    ASSERT_EQ(1, dist.size());
    ASSERT_NEAR(d * d, dist[0], 1e-6);
  }
}

