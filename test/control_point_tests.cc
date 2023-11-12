//
// (C) Matti Lehtonen 2023
//

#include "control_point.hpp"

#include <gtest/gtest.h>

#include <iostream>

using real = float;
using Point2 = curve::Point2<real>;
using Point3 = curve::Point3<real>;
using ControlPoint2 = curve::bezier::rational::ControlPoint<Point2>;
using ControlPoint3 = curve::bezier::rational::ControlPoint<Point3>;

struct ControlPointConstructorData {
  real const weight;
  Point3 const p;

  inline friend std::ostream &operator<<(std::ostream &out, ControlPointConstructorData const &data) {
    out << "{";
    out << data.weight << "," << data.p;
    out << "}";
    return out;
  }
};

class ControlPointConstructorTests : public ::testing::TestWithParam<ControlPointConstructorData> {};

TEST_P(ControlPointConstructorTests, Constructor2_correctlyConstructed) {
  auto const &params = GetParam();

  ControlPoint2 const result(params.weight, Point2(params.p.x(), params.p.y()));

  EXPECT_EQ(result.w(), params.weight);
  EXPECT_EQ(result.p().x(), params.p.x());
  EXPECT_EQ(result.p().y(), params.p.y());
}

TEST_P(ControlPointConstructorTests, Constructor3_correctlyConstructed) {
  auto const &params = GetParam();

  ControlPoint3 const result(params.weight, params.p);

  EXPECT_EQ(result.w(), params.weight);
  EXPECT_EQ(result.p().x(), params.p.x());
  EXPECT_EQ(result.p().y(), params.p.y());
  EXPECT_EQ(result.p().z(), params.p.z());
}

static constexpr ControlPointConstructorData controlPointConstructorData[]{
    {real(1), Point3(real(1.5), real(2.5), real(3.5))},
};

INSTANTIATE_TEST_SUITE_P(Fixture, ControlPointConstructorTests, ::testing::ValuesIn(controlPointConstructorData));
