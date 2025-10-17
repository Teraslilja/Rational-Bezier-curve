//
// (C) Matti Lehtonen 2023
//

#include "point2.hpp"
#include "point3.hpp"

#include <gtest/gtest.h>

#include <iostream>

using real = float;
using Point2 = curve::points::Point2<real>;
using Point3 = curve::points::Point3<real>;

TEST(PointConstructorTests, DefaultConstructor2_zeroValues)
{
  Point2 const result;

  EXPECT_EQ(result.x(), real(0));
  EXPECT_EQ(result.y(), real(0));
}

TEST(PointConstructorTests, DefaultConstructor3_zeroValues)
{
  Point3 const result;

  EXPECT_EQ(result.x(), real(0));
  EXPECT_EQ(result.y(), real(0));
  EXPECT_EQ(result.z(), real(0));
}

TEST(PointOrderTests, Order_2)
{
  EXPECT_EQ(Point2::getOrder(), 2u);
}

TEST(PointOrderTests, Order_3)
{
  EXPECT_EQ(Point3::getOrder(), 3u);
}

struct PointConstructorData
{
  real const x;
  real const y;
  real const z;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointConstructorData const& data)
  {
    out << "{";
    out << data.x << "," << data.y << "," << data.z;
    out << "}";
    return out;
  }
};

class PointConstructorTests
  : public ::testing::TestWithParam<PointConstructorData>
{};

TEST_P(PointConstructorTests, Constructor2_correctlyConstructed)
{
  auto const& params = GetParam();

  Point2 const result(params.x, params.y);

  EXPECT_EQ(result.x(), params.x);
  EXPECT_EQ(result.y(), params.y);
}

TEST_P(PointConstructorTests, Constructor3_correctlyConstructed)
{
  auto const& params = GetParam();

  Point3 const result(params.x, params.y, params.z);

  EXPECT_EQ(result.x(), params.x);
  EXPECT_EQ(result.y(), params.y);
  EXPECT_EQ(result.z(), params.z);
}

TEST_P(PointConstructorTests, Constructor3by2_correctlyConstructed)
{
  auto const& params = GetParam();

  Point3 const result(Point2(params.x, params.y));

  EXPECT_EQ(result.x(), params.x);
  EXPECT_EQ(result.y(), params.y);
  EXPECT_EQ(result.z(), real(0));
}

static constexpr PointConstructorData pointConstructorData[]{
  { real(1), real(2), real(3) },
  { real(-1), real(-2), real(-3) },
  { real(1.5), real(2.5), real(3.5) },
  { real(-1.5), real(-2.5), real(-3.5) },
};

INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointConstructorTests,
                         ::testing::ValuesIn(pointConstructorData));

struct PointOperatorMultiplyData
{
  Point3 const xyz;
  real const multi;

  Point2 const expectedResult2;
  Point3 const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointOperatorMultiplyData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.multi;
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

struct PointOperatorDivideData
{
  Point3 const xyz;
  real const divide;

  std::optional<Point2> const expectedResult2;
  std::optional<Point3> const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointOperatorDivideData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.divide << ",";
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

struct PointOperatorSubstractData
{
  Point3 const xyz;
  Point3 const substract;

  Point2 const expectedResult2;
  Point3 const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointOperatorSubstractData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.substract << ",";
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

struct PointOperatorAddData
{
  Point3 const xyz;
  Point3 const add;

  Point2 const expectedResult2;
  Point3 const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointOperatorAddData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.add << ",";
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

struct PointTraceData
{
  Point3 const xyz;

  real const expectedResult2;
  real const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointTraceData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

struct PointLengthSquaredData
{
  Point3 const xyz;

  real const expectedResult2;
  real const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointLengthSquaredData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

struct PointLengthData
{
  Point3 const xyz;

  real const expectedResult2;
  real const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointLengthData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

struct PointDistanceData
{
  Point3 const xyz;
  Point3 const p;

  real const expectedResult2;
  real const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointDistanceData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.p << ",";
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

struct PointNormalizedData
{
  Point3 const xyz;

  std::optional<Point2> const expectedResult2;
  std::optional<Point3> const expectedResult3;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         PointNormalizedData const& data)
  {
    out << "{";
    out << data.xyz << ",";
    out << data.expectedResult2 << "," << data.expectedResult3;
    out << "}";
    return out;
  }
};

class PointOperatorMultiplyTests
  : public ::testing::TestWithParam<PointOperatorMultiplyData>
{};
class PointOperatorDivideTests
  : public ::testing::TestWithParam<PointOperatorDivideData>
{};
class PointOperatorSubstractTests
  : public ::testing::TestWithParam<PointOperatorSubstractData>
{};
class PointOperatorAddTests
  : public ::testing::TestWithParam<PointOperatorAddData>
{};
class PointTraceTests : public ::testing::TestWithParam<PointTraceData>
{};
class PointLengthSquaredTests
  : public ::testing::TestWithParam<PointLengthSquaredData>
{};
class PointLengthTests : public ::testing::TestWithParam<PointLengthData>
{};
class PointDistanceTests : public ::testing::TestWithParam<PointDistanceData>
{};
class PointNormalizedTests
  : public ::testing::TestWithParam<PointNormalizedData>
{};

TEST_P(PointOperatorMultiplyTests, operatorMultiply2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp(params.xyz.x(), params.xyz.y());

  Point2 const result = tmp * params.multi;

  EXPECT_EQ(result.x(), params.expectedResult2.x());
  EXPECT_EQ(result.y(), params.expectedResult2.y());
}

TEST_P(PointOperatorMultiplyTests, operatorMultiply3_correctResult)
{
  auto const& params = GetParam();

  Point3 const result = params.xyz * params.multi;

  EXPECT_EQ(result.x(), params.expectedResult3.x());
  EXPECT_EQ(result.y(), params.expectedResult3.y());
  EXPECT_EQ(result.z(), params.expectedResult3.z());
}

TEST_P(PointOperatorDivideTests, operatorDivide2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp(params.xyz.x(), params.xyz.y());

  std::optional<Point2> const result = tmp / params.divide;

  ASSERT_EQ(result.has_value(), params.expectedResult2.has_value());
  if (params.expectedResult2.has_value()) {
    EXPECT_EQ(result.value().x(), params.expectedResult2.value().x());
    EXPECT_EQ(result.value().y(), params.expectedResult2.value().y());
  }
}

TEST_P(PointOperatorDivideTests, operatorDivide3_correctResult)
{
  auto const& params = GetParam();

  std::optional<Point3> const result = params.xyz / params.divide;

  ASSERT_EQ(result.has_value(), params.expectedResult3.has_value());
  if (params.expectedResult3.has_value()) {
    EXPECT_EQ(result.value().x(), params.expectedResult3.value().x());
    EXPECT_EQ(result.value().y(), params.expectedResult3.value().y());
    EXPECT_EQ(result.value().z(), params.expectedResult3.value().z());
  }
}

TEST_P(PointOperatorSubstractTests, operatorSubstract2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp1(params.xyz.x(), params.xyz.y());
  Point2 const tmp2(params.substract.x(), params.substract.y());

  Point2 const result = tmp1 - tmp2;

  EXPECT_EQ(result.x(), params.expectedResult2.x());
  EXPECT_EQ(result.y(), params.expectedResult2.y());
}

TEST_P(PointOperatorSubstractTests, operatorSubstract3_correctResult)
{
  auto const& params = GetParam();

  Point3 const result = params.xyz - params.substract;

  EXPECT_EQ(result.x(), params.expectedResult3.x());
  EXPECT_EQ(result.y(), params.expectedResult3.y());
  EXPECT_EQ(result.z(), params.expectedResult3.z());
}

TEST_P(PointOperatorAddTests, operatorAdd2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp2(params.add.x(), params.add.y());

  Point2 result(params.xyz.x(), params.xyz.y());
  result += tmp2;

  EXPECT_EQ(result.x(), params.expectedResult2.x());
  EXPECT_EQ(result.y(), params.expectedResult2.y());
}

TEST_P(PointOperatorAddTests, operatorAdd3_correctResult)
{
  auto const& params = GetParam();

  Point3 result(params.xyz.x(), params.xyz.y(), params.xyz.z());
  result += params.add;

  EXPECT_EQ(result.x(), params.expectedResult3.x());
  EXPECT_EQ(result.y(), params.expectedResult3.y());
  EXPECT_EQ(result.z(), params.expectedResult3.z());
}

TEST_P(PointOperatorAddTests, operatorSum2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp1(params.xyz.x(), params.xyz.y());
  Point2 const tmp2(params.add.x(), params.add.y());

  Point2 const result = tmp1 + tmp2;

  EXPECT_EQ(result.x(), params.expectedResult2.x());
  EXPECT_EQ(result.y(), params.expectedResult2.y());
}

TEST_P(PointOperatorAddTests, operatorSum3_correctResult)
{
  auto const& params = GetParam();

  Point3 const result = params.xyz + params.add;

  EXPECT_EQ(result.x(), params.expectedResult3.x());
  EXPECT_EQ(result.y(), params.expectedResult3.y());
  EXPECT_EQ(result.z(), params.expectedResult3.z());
}

TEST_P(PointTraceTests, trace2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp(params.xyz.x(), params.xyz.y());

  real const result = tmp.trace();

  EXPECT_EQ(result, params.expectedResult2);
}

TEST_P(PointTraceTests, trace3_correctResult)
{
  auto const& params = GetParam();

  real const result = params.xyz.trace();

  EXPECT_EQ(result, params.expectedResult3);
}

TEST_P(PointLengthSquaredTests, lengthSquared2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp(params.xyz.x(), params.xyz.y());

  real const result = tmp.lengthSquared();

  EXPECT_EQ(result, params.expectedResult2);
}

TEST_P(PointLengthSquaredTests, lengthSquared3_correctResult)
{
  auto const& params = GetParam();

  real const result = params.xyz.lengthSquared();

  EXPECT_EQ(result, params.expectedResult3);
}

TEST_P(PointLengthTests, length2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp(params.xyz.x(), params.xyz.y());

  real const result = tmp.length();

  EXPECT_EQ(result, params.expectedResult2);
}

TEST_P(PointLengthTests, length3_correctResult)
{
  auto const& params = GetParam();

  real const result = params.xyz.length();

  EXPECT_EQ(result, params.expectedResult3);
}

TEST_P(PointDistanceTests, distance2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp1(params.xyz.x(), params.xyz.y());
  Point2 const tmp2(params.p.x(), params.p.y());

  real const result = tmp1.distance(tmp2);

  EXPECT_EQ(result, params.expectedResult2);
}

TEST_P(PointDistanceTests, distance3_correctResult)
{
  auto const& params = GetParam();

  real const result = params.xyz.distance(params.p);

  EXPECT_EQ(result, params.expectedResult3);
}

TEST_P(PointNormalizedTests, normalized2_correctResult)
{
  auto const& params = GetParam();
  Point2 const tmp(params.xyz.x(), params.xyz.y());

  std::optional<Point2> const result = tmp.normalize();

  ASSERT_EQ(result.has_value(), params.expectedResult2.has_value());
  if (params.expectedResult2.has_value()) {
    EXPECT_EQ(result.value().x(), params.expectedResult2.value().x());
    EXPECT_EQ(result.value().y(), params.expectedResult2.value().y());
  }
}

TEST_P(PointNormalizedTests, normalized3_correctResult)
{
  auto const& params = GetParam();

  std::optional<Point3> const result = params.xyz.normalize();

  ASSERT_EQ(result.has_value(), params.expectedResult3.has_value());
  if (params.expectedResult3.has_value()) {
    EXPECT_EQ(result.value().x(), params.expectedResult3.value().x());
    EXPECT_EQ(result.value().y(), params.expectedResult3.value().y());
    EXPECT_EQ(result.value().z(), params.expectedResult3.value().z());
  }
}

static constexpr PointOperatorMultiplyData pointOperatorMultiplyData[]{
  { Point3(real(1.5), real(2.5), real(3.5)),
    real(2),
    Point2(real(3), real(5)),
    Point3(real(3), real(5), real(7)) },
};

static constexpr PointOperatorDivideData pointOperatorDivideData[]{
  { Point3(real(2), real(3), real(4)),
    real(2),
    Point2(real(1), real(1.5)),
    Point3(real(1), real(1.5), real(2)) },
  { Point3(real(2), real(3), real(4)), real(0), std::nullopt, std::nullopt },
  { Point3(real(0), real(0), real(0)), real(0), std::nullopt, std::nullopt },
};

static constexpr PointOperatorSubstractData pointOperatorSubstractData[]{
  { Point3(real(4), real(6), real(8)),
    Point3(real(1), real(2), real(3)),
    Point2(real(3), real(4)),
    Point3(real(3), real(4), real(5)) },
};

static constexpr PointOperatorAddData pointOperatorAddData[]{
  { Point3(real(4), real(6), real(8)),
    Point3(real(1), real(2), real(3)),
    Point2(real(5), real(8)),
    Point3(real(5), real(8), real(11)) },
};

static constexpr PointTraceData pointTraceData[]{
  { Point3(real(-4), real(6), real(1)), real(2), real(3) },
  { Point3(real(4), real(-6), real(1)), real(-2), real(-1) },
  { Point3(real(4), real(6), real(-1)), real(10), real(9) },
};

static constexpr PointLengthSquaredData pointLengthSquaredData[]{
  { Point3(real(3), real(4), real(5)), real(25), real(50) },
};

static const PointLengthData pointLengthData[]{
  { Point3(real(3), real(4), real(5)), real(5), real(std::sqrt(50)) },
};

static constexpr PointDistanceData pointDistanceData[]{
  { Point3(real(1), real(1), real(1)),
    Point3(real(4), real(5), real(1)),
    real(5),
    real(5) },
};

static inline constexpr real one_over_sqrt_2 =
  real(1) / std::numbers::sqrt2_v<real>;
static inline constexpr real one_over_sqrt_3 = std::numbers::inv_sqrt3_v<real>;

static constexpr PointNormalizedData pointNormalizedData[]{
  { Point3(real(0), real(0), real(0)), std::nullopt, std::nullopt },
  { Point3(real(0), real(0), real(2)),
    std::nullopt,
    Point3(real(0), real(0), real(1)) },
  { Point3(real(1), real(1), real(1)),
    Point2(one_over_sqrt_2, one_over_sqrt_2),

    Point3(one_over_sqrt_3, one_over_sqrt_3, one_over_sqrt_3) },
};

INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointOperatorMultiplyTests,
                         ::testing::ValuesIn(pointOperatorMultiplyData));
INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointOperatorDivideTests,
                         ::testing::ValuesIn(pointOperatorDivideData));
INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointOperatorSubstractTests,
                         ::testing::ValuesIn(pointOperatorSubstractData));
INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointOperatorAddTests,
                         ::testing::ValuesIn(pointOperatorAddData));
INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointTraceTests,
                         ::testing::ValuesIn(pointTraceData));
INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointLengthSquaredTests,
                         ::testing::ValuesIn(pointLengthSquaredData));
INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointLengthTests,
                         ::testing::ValuesIn(pointLengthData));
INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointDistanceTests,
                         ::testing::ValuesIn(pointDistanceData));
INSTANTIATE_TEST_SUITE_P(Fixture,
                         PointNormalizedTests,
                         ::testing::ValuesIn(pointNormalizedData));
