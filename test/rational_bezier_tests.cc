//
// (C) Matti Lehtonen 2023
//

#include "control_point.hpp"
#include "rational_bezier.hpp"

#include <gtest/gtest.h>

#include <iostream>

using real = float;
using Point2 = curve::Point2<real>;
using Point3 = curve::Point3<real>;
using ControlPoint2 = curve::bezier::rational::ControlPoint<Point2>;
using ControlPoint3 = curve::bezier::rational::ControlPoint<Point3>;
using Rational2 = curve::bezier::rational::Rational<ControlPoint2>;
using curve::bezier::rational::ValidityIssue;

using real_or_issue = Rational2::real_or_issue;

using Point2_or_issue = Rational2::point_or_issue;
using vector_of_Point2s_or_issue = Rational2::vector_of_points_or_issue;

static constexpr ControlPoint2 cp100 = {real(1), Point2(real(0), real(0))};
static constexpr ControlPoint2 cp110 = {real(1), Point2(real(1), real(0))};
static constexpr ControlPoint2 cp101 = {real(1), Point2(real(0), real(1))};
static constexpr ControlPoint2 cp111 = {real(1), Point2(real(1), real(1))};

static constexpr std::vector<ControlPoint2> empty{};
static const std::vector<ControlPoint2> single_point{cp100};

/* | */
static const std::vector<ControlPoint2> vertical_line{cp100, cp101};
/* _ */
static const std::vector<ControlPoint2> horizontal_line{cp100, cp110};
/* / */
static const std::vector<ControlPoint2> slope_line{cp100, cp111};

static constexpr real triangle_side = real(1);
static const real triangle_height = std::sqrt(real(0.75) * triangle_side);
static constexpr ControlPoint2 cp_arc0 = {real(1), Point2(real(0) * triangle_side, real(0))};
static const ControlPoint2 cp_arc1 = {real(0.5), Point2(real(0.5) * triangle_side, triangle_height)};
static constexpr ControlPoint2 cp_arc2end = {real(1), Point2(real(1) * triangle_side, real(0))};
static const ControlPoint2 cp_arc2 = {real(0.5), Point2(real(1) * triangle_side, real(0))};
static const ControlPoint2 cp_arc3 = {real(0.5), Point2(real(1.5) * triangle_side, triangle_height)};
static constexpr ControlPoint2 cp_arc4end = {real(1), Point2(real(2) * triangle_side, real(0))};

// Curve with bad weights
static constexpr ControlPoint2 cp000 = {real(0), Point2(real(0), real(0))};
static constexpr ControlPoint2 cp011 = {real(0), Point2(real(1), real(1))};
static const std::vector<ControlPoint2> bad_curve{cp000, cp011};

/* /\ */
static const std::vector<ControlPoint2> single_arc{cp_arc0, cp_arc1, cp_arc2end};
/* /\/\ */
static const std::vector<ControlPoint2> three_arcs{cp_arc0, cp_arc1, cp_arc2, cp_arc3, cp_arc4end};

template <class T> inline std::ostream &operator<<(std::ostream &out, std::variant<T, ValidityIssue> const &data) {
  if (std::holds_alternative<T>(data)) {
    out << std::get<T>(data);
  }
  if (std::holds_alternative<ValidityIssue>(data)) {
    out << std::get<ValidityIssue>(data);
  }
  return out;
}

template <class T> inline std::ostream &operator<<(std::ostream &out, std::optional<T> const &data) {
  if (data.has_value()) {
    out << data.value();
  } else {
    out << "'no value'";
  }
  return out;
}

// Helpers
static void validateResult(real_or_issue const &result, real_or_issue const &expected, const real EPS) {
  ASSERT_EQ(std::holds_alternative<real>(result), std::holds_alternative<real>(expected));
  ASSERT_EQ(std::holds_alternative<ValidityIssue>(result), std::holds_alternative<ValidityIssue>(expected));
  if (std::holds_alternative<real>(expected)) {
    EXPECT_NEAR(std::get<real>(result), std::get<real>(expected), EPS);
  }
  if (std::holds_alternative<ValidityIssue>(expected)) {
    EXPECT_EQ(std::get<ValidityIssue>(result), std::get<ValidityIssue>(expected));
  }
}

static void validateResult(Point2_or_issue const &result, Point2_or_issue const &expected, const real EPS) {
  ASSERT_EQ(std::holds_alternative<Point2>(result), std::holds_alternative<Point2>(expected));
  ASSERT_EQ(std::holds_alternative<ValidityIssue>(result), std::holds_alternative<ValidityIssue>(expected));
  if (std::holds_alternative<Point2>(expected)) {
    EXPECT_NEAR(std::get<Point2>(result).x(), std::get<Point2>(expected).x(), EPS);
    EXPECT_NEAR(std::get<Point2>(result).y(), std::get<Point2>(expected).y(), EPS);
  }
  if (std::holds_alternative<ValidityIssue>(expected)) {
    EXPECT_EQ(std::get<ValidityIssue>(result), std::get<ValidityIssue>(expected));
  }
}

struct RationalConstructorData {
  std::vector<ControlPoint2> const cp;

  inline friend std::ostream &operator<<(std::ostream &out, RationalConstructorData const &data) {
    out << "{";
    out << data.cp;
    out << "}";
    return out;
  }
};

class RationalConstructorTests : public ::testing::TestWithParam<RationalConstructorData> {};

TEST_P(RationalConstructorTests, Constructor_correctlyConstructed) {
  auto const &params = GetParam();

  Rational2 const result{params.cp};

  ASSERT_EQ(result.numberOfControlPoints(), params.cp.size());
  for (std::size_t i = 0u; i < params.cp.size(); ++i) {
    auto const has_cp = result.getControlPoint(i);
    ASSERT_TRUE(has_cp.has_value());
    ControlPoint2 const &cp = has_cp.value();
    EXPECT_EQ(cp.w(), params.cp.at(i).w());
    EXPECT_EQ(cp.p().x(), params.cp.at(i).p().x());
    EXPECT_EQ(cp.p().y(), params.cp.at(i).p().y());
  }
}

static const RationalConstructorData rationalConstructorData[]{
    {empty},      {single_point}, {vertical_line}, {horizontal_line},
    {slope_line}, {single_arc},   {three_arcs},    {bad_curve},
};

INSTANTIATE_TEST_SUITE_P(Fixture, RationalConstructorTests, ::testing::ValuesIn(rationalConstructorData));

struct PointAtData {
  Rational2 const curve;
  real const u;

  Point2_or_issue const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, PointAtData const &data) {
    out << "{";
    out << data.curve << "," << data.u << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class PointAtTests : public ::testing::TestWithParam<PointAtData> {};

TEST_P(PointAtTests, pointAt_correctValue) {
  static constexpr real EPS = real(1e-5);
  auto const &params = GetParam();

  Point2_or_issue const result = params.curve.pointAt(params.u);

  validateResult(result, params.expectedResult, EPS);
}

static const PointAtData pointAtData[]{
    {Rational2(empty), real(0), ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS},
    {Rational2(single_point), real(-1), ValidityIssue::ISSUE_U_IS_INVALID},
    {Rational2(single_point), real(2), ValidityIssue::ISSUE_U_IS_INVALID},

    {Rational2(single_point), real(0), cp100.p()},
    {Rational2(single_point), real(1), cp100.p()},
    {Rational2(vertical_line), real(0), cp100.p()},
    {Rational2(vertical_line), real(0.5), ((cp100.p() + cp101.p()) / real(2)).value()},
    {Rational2(vertical_line), real(1), cp101.p()},
    {Rational2(horizontal_line), real(0), cp100.p()},
    {Rational2(horizontal_line), real(0.5), ((cp100.p() + cp110.p()) / real(2)).value()},
    {Rational2(horizontal_line), real(1), cp110.p()},
    {Rational2(slope_line), real(0), cp100.p()},
    {Rational2(slope_line), real(0.5), ((cp100.p() + cp111.p()) / real(2)).value()},
    {Rational2(slope_line), real(1), cp111.p()},
    {Rational2(single_arc), real(0), cp_arc0.p()},
    {Rational2(single_arc), real(0.5), Point2(real(0.5) * triangle_side, real(0.288675))},
    {Rational2(single_arc), real(1), cp_arc2end.p()},
    {Rational2(three_arcs), real(0), cp_arc0.p()},
    {Rational2(three_arcs), real(0.5), Point2(triangle_side, real(0.384900))},
    {Rational2(three_arcs), real(1), cp_arc4end.p()},
};

INSTANTIATE_TEST_SUITE_P(Fixture, PointAtTests, ::testing::ValuesIn(pointAtData));

struct CurveLengthData {
  Rational2 const curve;

  real_or_issue const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, CurveLengthData const &data) {
    out << "{";
    out << data.curve << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class CurveLengthTests : public ::testing::TestWithParam<CurveLengthData> {};

TEST_P(CurveLengthTests, curveLength_correctValue) {
  static constexpr real EPS = real(1e-5);
  auto const &params = GetParam();

  real_or_issue const result = params.curve.curveLength();

  validateResult(result, params.expectedResult, EPS);
}

static const CurveLengthData curveLengthData[]{
    {Rational2(empty), ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS},
    {Rational2(single_point), real(0)},
    {Rational2(vertical_line), real(1)},
    {Rational2(horizontal_line), real(1)},
    {Rational2(slope_line), std::sqrt(real(2))},
    {Rational2(single_arc), real(1.209144)},
    {Rational2(three_arcs), real(2.253819)},

    {Rational2(bad_curve), ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS},
};

INSTANTIATE_TEST_SUITE_P(Fixture, CurveLengthTests, ::testing::ValuesIn(curveLengthData));

struct VelocitySpeedTangentData {
  Rational2 const curve;
  real const u;

  Point2_or_issue const expectedVelocityResult;
  real_or_issue const expectedSpeedResult;
  Point2_or_issue const expectedTangentResult;

  inline friend std::ostream &operator<<(std::ostream &out, VelocitySpeedTangentData const &data) {
    out << "{";
    out << data.curve << "," << data.u << ",";
    out << data.expectedVelocityResult << "," << data.expectedSpeedResult << "," << data.expectedTangentResult;
    out << "}";
    return out;
  }
};

class VelocitySpeedTangentTests : public ::testing::TestWithParam<VelocitySpeedTangentData> {
  [[nodiscard]] static constexpr std::optional<Point2> approximateTangentAt(Rational2 const &curve,
                                                                            real const u) noexcept {
    if ((u < real(0)) || (u > real(1)))
      return std::nullopt;

    switch (curve.numberOfControlPoints()) {
    case 0u:
    case 1u:
      return std::nullopt;

    default:
      // Approximate tangent
      constexpr real du = real(1e-5);
      Point2_or_issue const has_p0 = curve.pointAt(std::max(real(0), u - du));
      if (!std::holds_alternative<Point2>(has_p0)) {
        return std::nullopt;
      }
      Point2_or_issue const has_p1 = curve.pointAt(std::min(real(1), u + du));
      if (!std::holds_alternative<Point2>(has_p1)) {
        return std::nullopt;
      }
      Point2 const delta = std::get<Point2>(has_p1) - std::get<Point2>(has_p0);
      return delta.normalize();
    }
  }
};

TEST_P(VelocitySpeedTangentTests, velocityAt_correctValue) {
  static constexpr real EPS = real(1e-5);
  auto const &params = GetParam();

  Point2_or_issue const result = params.curve.velocityAt(params.u);

  validateResult(result, params.expectedVelocityResult, EPS);
}

TEST_P(VelocitySpeedTangentTests, speedAt_correctValue) {
  static constexpr real EPS = real(1e-5);
  auto const &params = GetParam();

  real_or_issue const result = params.curve.speedAt(params.u);

  validateResult(result, params.expectedSpeedResult, EPS);
}

TEST_P(VelocitySpeedTangentTests, tangentAt_correctValue) {
  static constexpr real EPS = real(1e-5);
  auto const &params = GetParam();

  Point2_or_issue const result = params.curve.tangentAt(params.u);

  validateResult(result, params.expectedTangentResult, EPS);
}

static const real one_per_sqrt2 = real(1) / std::sqrt(real(2));

static const VelocitySpeedTangentData velocitySpeedTangentData[]{
    {Rational2(empty), real(0.5), ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS,
     ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS, ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS},
    {Rational2(single_point), real(0.5), ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS,
     ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS, ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS},
    {Rational2(vertical_line), real(-1), ValidityIssue::ISSUE_U_IS_INVALID, ValidityIssue::ISSUE_U_IS_INVALID,
     ValidityIssue::ISSUE_U_IS_INVALID},
    {Rational2(vertical_line), real(2), ValidityIssue::ISSUE_U_IS_INVALID, ValidityIssue::ISSUE_U_IS_INVALID,
     ValidityIssue::ISSUE_U_IS_INVALID},

    {Rational2(bad_curve), real(0.5), ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS,
     ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS, ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS},

    {Rational2(vertical_line), real(0), Point2(real(0), real(1)), real(1), Point2(real(0), real(1))},
    {Rational2(vertical_line), real(0.5), Point2(real(0), real(1)), real(1), Point2(real(0), real(1))},
    {Rational2(vertical_line), real(1), Point2(real(0), real(1)), real(1), Point2(real(0), real(1))},

    {Rational2(horizontal_line), real(0), Point2(real(1), real(0)), real(1), Point2(real(1), real(0))},
    {Rational2(horizontal_line), real(0.5), Point2(real(1), real(0)), real(1), Point2(real(1), real(0))},
    {Rational2(horizontal_line), real(1), Point2(real(1), real(0)), real(1), Point2(real(1), real(0))},

    {Rational2(slope_line), real(0), Point2(real(1), real(1)), std::sqrt(real(2)),
     Point2(one_per_sqrt2, one_per_sqrt2)},
    {Rational2(slope_line), real(0.5), Point2(real(1), real(1)), std::sqrt(real(2)),
     Point2(one_per_sqrt2, one_per_sqrt2)},
    {Rational2(slope_line), real(1), Point2(real(1), real(1)), std::sqrt(real(2)),
     Point2(one_per_sqrt2, one_per_sqrt2)},

    {Rational2(single_arc), real(0), Point2(real(0.5), real(0.866025)), real(1), Point2(real(0.5), real(0.866025))},
    {Rational2(single_arc), real(0.5), Point2(real(1.3333331), real(0)), real(1.333333), Point2(real(1), real(0))},
    {Rational2(single_arc), real(1), Point2(real(0.5), real(-0.866025)), real(1),
     Point2(real(0.5), real(-0.866025))},

    //{Rational2(three_arcs),std::make_optional(real(2.253819))},
};

INSTANTIATE_TEST_SUITE_P(Fixture, VelocitySpeedTangentTests, ::testing::ValuesIn(velocitySpeedTangentData));

struct ControlPointCountData {
  Rational2 const curve;

  std::size_t const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, ControlPointCountData const &data) {
    out << "{";
    out << data.curve << ",";
    out << data.expectedResult;
    out << "}";
    return out;
  }
};

class ControlPointCountTests : public ::testing::TestWithParam<ControlPointCountData> {};

TEST_P(ControlPointCountTests, numberOfControlPoints_correctValue) {
  auto const &params = GetParam();

  std::size_t const result = params.curve.numberOfControlPoints();

  EXPECT_EQ(result, params.expectedResult);
}

static const ControlPointCountData controlPointCountData[]{
    {Rational2(empty), 0u},           {Rational2(single_point), 1u}, {Rational2(vertical_line), 2u},
    {Rational2(horizontal_line), 2u}, {Rational2(slope_line), 2u},   {Rational2(single_arc), 3u},
    {Rational2(three_arcs), 5u},
};

INSTANTIATE_TEST_SUITE_P(Fixture, ControlPointCountTests, ::testing::ValuesIn(controlPointCountData));

struct GetModifyAddRemoveControlPointData {
  std::vector<ControlPoint2> const controlPoints;
  std::size_t const index;

  std::optional<ControlPoint2> const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, GetModifyAddRemoveControlPointData const &data) {
    out << "{";
    out << data.controlPoints << "," << data.index << ",";
    out << data.expectedResult;
    out << "}";
    return out;
  }
};

class GetModifyAddRemoveControlPointTests : public ::testing::TestWithParam<GetModifyAddRemoveControlPointData> {};

TEST_P(GetModifyAddRemoveControlPointTests, getControlPoint_correctValue) {
  auto const &params = GetParam();
  Rational2 const curve(params.controlPoints);

  auto const result = curve.getControlPoint(params.index);

  ASSERT_EQ(result.has_value(), params.expectedResult.has_value());
  if (params.expectedResult.has_value()) {
    EXPECT_EQ(result.value().w(), params.expectedResult.value().w());
    EXPECT_EQ(result.value().p().x(), params.expectedResult.value().p().x());
    EXPECT_EQ(result.value().p().y(), params.expectedResult.value().p().y());
  }
}

TEST_P(GetModifyAddRemoveControlPointTests, modifyControlPoint_correctValue) {
  auto const &params = GetParam();
  Rational2 curve(params.controlPoints);
  auto cp = curve.getControlPoint(params.index);
  if (cp.has_value()) {
    real constexpr RESET_WEIGHT = real(0);
    real constexpr RESET_X = real(-1);
    real constexpr RESET_Y = real(-2);
    ASSERT_NE(cp.value().get().w(), RESET_WEIGHT);
    ASSERT_NE(cp.value().get().p().x(), RESET_X);
    ASSERT_NE(cp.value().get().p().y(), RESET_Y);

    cp.value().get().w() = RESET_WEIGHT;
    cp.value().get().p().x() = RESET_X;
    cp.value().get().p().y() = RESET_Y;
    auto const result = curve.getControlPoint(params.index);

    EXPECT_EQ(result.value().get().w(), RESET_WEIGHT);
    EXPECT_EQ(result.value().get().p().x(), RESET_X);
    EXPECT_EQ(result.value().get().p().y(), RESET_Y);
  }
}

TEST_P(GetModifyAddRemoveControlPointTests, removeControlPoint_correctValue) {
  auto const &params = GetParam();
  Rational2 curve(params.controlPoints);
  std::size_t const oldCount = curve.numberOfControlPoints();
  auto const cp = curve.getControlPoint(params.index);
  if (cp.has_value()) {
    ControlPoint2 const removed = cp.value().get();

    bool const result = curve.removeControlPoint(params.index);
    std::size_t const newCount = curve.numberOfControlPoints();

    EXPECT_EQ(newCount + 1u, oldCount);
    EXPECT_TRUE(result);
    auto const cp2 = curve.getControlPoint(params.index);
    if (cp2.has_value()) {
      EXPECT_NE(cp2.value().get().p().x(), removed.p().x());
      EXPECT_NE(cp2.value().get().p().y(), removed.p().y());
    }
  } else {

    bool const result = curve.removeControlPoint(params.index);
    std::size_t const newCount = curve.numberOfControlPoints();

    EXPECT_EQ(newCount, oldCount);
    EXPECT_FALSE(result);
  }
}

TEST_P(GetModifyAddRemoveControlPointTests, addControlPoint_correctValue) {
  static real constexpr RESET_WEIGHT = real(0);
  static real constexpr RESET_X = real(-1);
  static real constexpr RESET_Y = real(-2);
  static constexpr ControlPoint2 addedCp(real(RESET_WEIGHT), Point2(RESET_X, RESET_Y));
  auto const &params = GetParam();
  Rational2 curve(params.controlPoints);
  std::size_t const oldCount = curve.numberOfControlPoints();
  auto const cp = curve.getControlPoint(params.index);
  std::optional<ControlPoint2> const moved = cp.has_value() ? std::make_optional(cp.value().get()) : std::nullopt;

  bool const result = curve.addControlPoint(params.index, addedCp);
  std::size_t const newCount = curve.numberOfControlPoints();

  EXPECT_EQ(newCount, oldCount + 1u);
  EXPECT_TRUE(result);
  auto const cp1 = curve.getControlPoint(params.index);
  ASSERT_TRUE(cp1.has_value());
  EXPECT_EQ(cp1.value().get().w(), RESET_WEIGHT);
  EXPECT_EQ(cp1.value().get().p().x(), RESET_X);
  EXPECT_EQ(cp1.value().get().p().y(), RESET_Y);
  auto const cp2 = curve.getControlPoint(params.index + 1u);
  ASSERT_EQ(cp.has_value(), cp2.has_value());
  if (cp2.has_value()) {
    EXPECT_EQ(cp2.value().get().w(), moved.value().w());
    EXPECT_EQ(cp2.value().get().p().x(), moved.value().p().x());
    EXPECT_EQ(cp2.value().get().p().y(), moved.value().p().y());
  }
}

TEST_P(GetModifyAddRemoveControlPointTests, removeAllControlPoints_correctValue) {
  auto const &params = GetParam();
  Rational2 curve(params.controlPoints);

  curve.removeAllControlPoints();
  std::size_t const result = curve.numberOfControlPoints();

  EXPECT_EQ(result, 0u);
}

TEST_P(GetModifyAddRemoveControlPointTests, releaseMemory_correctValue) {
  auto const &params = GetParam();
  Rational2 curve(params.controlPoints);
  curve.removeAllControlPoints();

  curve.slim();
  std::size_t const result = curve.capacity();

  EXPECT_EQ(result, 0u);
}

TEST_P(GetModifyAddRemoveControlPointTests, reserveMemory_correctValue) {
  std::size_t constexpr ALOT = 1000u;
  auto const &params = GetParam();
  Rational2 curve(params.controlPoints);
  ASSERT_NE(curve.capacity(), ALOT);

  bool const result = curve.reserve(ALOT);

  ASSERT_TRUE(result); // Allocation is always succesful, no stress testing
  EXPECT_EQ(curve.capacity(), ALOT);
}

TEST_P(GetModifyAddRemoveControlPointTests, reserveMemory_fails) {
  std::size_t constexpr IMPOSSIBLE_AMOUNT = -1;
  auto const &params = GetParam();
  Rational2 curve(params.controlPoints);

  bool const result = curve.reserve(IMPOSSIBLE_AMOUNT);

  ASSERT_FALSE(result); // Allocation is always failure
}

static const GetModifyAddRemoveControlPointData getModifyAddRemoveControlPointData[]{
    {empty, 0u, std::nullopt}, {single_point, 0u, cp100},    {single_point, 1u, std::nullopt},
    {three_arcs, 0u, cp_arc0}, {three_arcs, 1u, cp_arc1},    {three_arcs, 2u, cp_arc2},
    {three_arcs, 3u, cp_arc3}, {three_arcs, 4u, cp_arc4end}, {three_arcs, 5u, std::nullopt},
};

INSTANTIATE_TEST_SUITE_P(Fixture, GetModifyAddRemoveControlPointTests,
                         ::testing::ValuesIn(getModifyAddRemoveControlPointData));

struct ClosestCurvePointData {
  Rational2 const curve;
  Point2 const point;

  Point2_or_issue const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, ClosestCurvePointData const &data) {
    out << "{";
    out << data.curve << "," << data.point << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class ClosestCurvePointTests : public ::testing::TestWithParam<ClosestCurvePointData> {};

TEST_P(ClosestCurvePointTests, getClosestPoint_correctValue) {
  static real constexpr EPS = real(1e-4);
  auto const &params = GetParam();

  auto const result = params.curve.closestCurvePointFor(params.point);

  validateResult(result, params.expectedResult, EPS);
}

static const ClosestCurvePointData closestCurvePointData[]{
    {Rational2(empty), Point2(real(0), real(0)), ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS},

    {Rational2(single_point), cp100.p(), cp100.p()},
    {Rational2(single_point), cp111.p(), cp100.p()},

    {Rational2(slope_line), cp100.p(), cp100.p()},
    {Rational2(slope_line), cp111.p(), cp111.p()},
    {Rational2(slope_line), cp110.p(), Point2(real(0.5), real(0.5))},
    {Rational2(slope_line), cp101.p(), Point2(real(0.5), real(0.5))},
    {Rational2(slope_line), Point2(real(0.5), real(0.5)), Point2(real(0.5), real(0.5))},

    {Rational2(single_arc), cp_arc1.p() + Point2(real(0), real(1)), Point2(real(0.5), real(0.288675))},

    {Rational2(bad_curve), Point2(real(0), real(0)), ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS},
};

INSTANTIATE_TEST_SUITE_P(Fixture, ClosestCurvePointTests, ::testing::ValuesIn(closestCurvePointData));

struct LineStringData {
  Rational2 const curve;

  std::optional<ValidityIssue> const expectedIssue;
  std::optional<std::size_t> const expectedSize;

  inline friend std::ostream &operator<<(std::ostream &out, LineStringData const &data) {
    out << "{";
    out << data.curve << "," << data.expectedIssue << "," << data.expectedSize;
    out << "}";
    return out;
  }
};

class LineStringTests : public ::testing::TestWithParam<LineStringData> {};

TEST_P(LineStringTests, asLineString_correctValue) {
  auto const &params = GetParam();

  vector_of_Point2s_or_issue const result = params.curve.asLinestring();

  ASSERT_EQ(std::holds_alternative<std::vector<Point2>>(result), params.expectedSize.has_value());
  ASSERT_EQ(std::holds_alternative<ValidityIssue>(result), params.expectedIssue.has_value());
  if (params.expectedSize.has_value()) {
    EXPECT_EQ(std::get<std::vector<Point2>>(result).size(), params.expectedSize.value());
  }
  if (params.expectedIssue.has_value()) {
    EXPECT_EQ(std::get<ValidityIssue>(result), params.expectedIssue.value());
  }
}

static const LineStringData lineStringData[]{
    {Rational2(empty), ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS, std::nullopt},
    {Rational2(single_point), std::nullopt, 1u},
    {Rational2(vertical_line), std::nullopt, 2u},
    {Rational2(horizontal_line), std::nullopt, 2u},
    {Rational2(slope_line), std::nullopt, 2u},
    {Rational2(single_arc), std::nullopt, 65u},
    {Rational2(three_arcs), std::nullopt, 101u},

    {Rational2(bad_curve), ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS, std::nullopt},
};

INSTANTIATE_TEST_SUITE_P(Fixture, LineStringTests, ::testing::ValuesIn(lineStringData));
