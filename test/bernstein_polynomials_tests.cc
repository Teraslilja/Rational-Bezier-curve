//
// (C) Matti Lehtonen 2023
//

#include "bernstein_polynomials.hpp"

#include <gtest/gtest.h>

#include <iostream>

using real = float;

struct Testing : public curve::bezier::utilities::BernsteinPolynomials<real> {
  using Base = curve::bezier::utilities::BernsteinPolynomials<real>;

  [[nodiscard]] static inline constexpr std::size_t binomial(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomial(n, k);
  }

  [[nodiscard]] static inline constexpr real toIntegerPower(real x, std::size_t i) noexcept {
    return Base::toIntegerPower(x, i);
  }
};

struct BinomialData {
  std::size_t const n;
  std::size_t const k;

  std::size_t const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, BinomialData const &data) {
    out << "{";
    out << data.n << "," << data.k << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class BinomialTests : public ::testing::TestWithParam<BinomialData> {};

TEST_P(BinomialTests, binomial_correctResult) {
  auto const &params = GetParam();

  std::size_t const result = Testing::binomial(params.n, params.k);

  EXPECT_EQ(result, params.expectedResult);
}

//      1
//     1 1
//    1 2   1
//   1 3   3 1
//  1 4  6  4 1
// 1 5 10 10 5 1
static constexpr BinomialData binomialData[]{
    {0u, 1u, 0u}, {0u, 0u, 1u},

    {1u, 0u, 1u}, {1u, 1u, 1u},

    {2u, 0u, 1u}, {2u, 1u, 2u}, {2u, 2u, 1u},

    {3u, 0u, 1u}, {3u, 1u, 3u}, {3u, 2u, 3u},  {3u, 3u, 1u},

    {4u, 0u, 1u}, {4u, 1u, 4u}, {4u, 2u, 6u},  {4u, 3u, 4u},  {4u, 4u, 1u},

    {5u, 0u, 1u}, {5u, 1u, 5u}, {5u, 2u, 10u}, {5u, 3u, 10u}, {5u, 4u, 5u}, {5u, 5u, 1u},
};

INSTANTIATE_TEST_SUITE_P(Fixture, BinomialTests, ::testing::ValuesIn(binomialData));

struct IntegerPowerData {
  real const x;
  std::size_t const i;

  real const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, IntegerPowerData const &data) {
    out << "{";
    out << data.x << "," << data.i << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class IntegerPowerTests : public ::testing::TestWithParam<IntegerPowerData> {};

TEST_P(IntegerPowerTests, integerPower_correctResult) {
  auto const &params = GetParam();

  real const result = Testing::toIntegerPower(params.x, params.i);

  EXPECT_EQ(result, params.expectedResult);
}

static constexpr IntegerPowerData integerPowerData[]{
    {real(0), 0u, real(1)}, {real(0.5), 0u, real(1)},     {real(1), 0u, real(1)},

    {real(0), 1u, real(0)}, {real(0.5), 1u, real(0.5)},   {real(1), 1u, real(1)},

    {real(0), 2u, real(0)}, {real(0.5), 2u, real(0.25)},  {real(1), 2u, real(1)},

    {real(0), 3u, real(0)}, {real(0.5), 3u, real(0.125)}, {real(1), 3u, real(1)},
};

INSTANTIATE_TEST_SUITE_P(Fixture, IntegerPowerTests, ::testing::ValuesIn(integerPowerData));

struct BernsteinPolynomialData {
  std::size_t const i;
  std::size_t const n;
  real const u;

  real const expectedResultB;
  real const expectedResultdB;

  inline friend std::ostream &operator<<(std::ostream &out, BernsteinPolynomialData const &data) {
    out << "{";
    out << data.i << "," << data.n << "," << data.u << "," << data.expectedResultB << "," << data.expectedResultdB;
    out << "}";
    return out;
  }
};

class BernsteinPolynomialTests : public ::testing::TestWithParam<BernsteinPolynomialData> {};

TEST_P(BernsteinPolynomialTests, B_correctResult) {
  auto const &params = GetParam();

  real const result = Testing::B(params.i, params.n, params.u);

  EXPECT_EQ(result, params.expectedResultB);
}

TEST_P(BernsteinPolynomialTests, dB_correctResult) {
  auto const &params = GetParam();

  real const result = Testing::dB(params.i, params.n, params.u);

  EXPECT_EQ(result, params.expectedResultdB);
}

static constexpr BernsteinPolynomialData bernsteinPolynomialData[]{
    {0u, 1u, real(0), real(1 * 1.0 * 1.0), real(1 * (0.0 - 1.0))},
    {0u, 1u, real(0.5), real(1 * 1.0 * 0.5), real(1 * (0.0 - 1.0))},
    {0u, 1u, real(1), real(1 * 1.0 * 0.0), real(1 * (0.0 - 1.0))},

    {0u, 2u, real(0), real(1 * 1.0 * 1.0), real(2 * (0.0 - 1.0))},
    {0u, 2u, real(0.5), real(1 * 1.0 * 0.25), real(2 * (0.0 - 0.5))},
    {0u, 2u, real(1), real(1 * 1.0 * 0.0), real(2 * (0.0 - 0.0))},

    {1u, 2u, real(0), real(2 * 0.0 * 1.0), real(2 * (1.0 - 0.0))},
    {1u, 2u, real(0.5), real(2 * 0.5 * 0.5), real(2 * (0.5 - 0.5))},
    {1u, 2u, real(1), real(2 * 1.0 * 0.0), real(2 * (0.0 - 1.0))},

    {1u, 3u, real(0), real(3 * 0.0 * 1.0), real(3 * (1.0 - 0.0))},
    {1u, 3u, real(0.5), real(3 * 0.5 * 0.25), real(3 * (0.25 - 0.5))},
    {1u, 3u, real(1), real(3 * 1.0 * 0.0), real(3 * (0.0 - 0.0))},

    {2u, 3u, real(0), real(3 * 0.0 * 1.0), real(3 * (0.0 - 0.0))},
    {2u, 3u, real(0.5), real(3 * 0.25 * 0.5), real(3 * (0.5 - 0.25))},
    {2u, 3u, real(1), real(3 * 1.0 * 0.0), real(3 * (0.0 - 1.0))},
};

INSTANTIATE_TEST_SUITE_P(Fixture, BernsteinPolynomialTests, ::testing::ValuesIn(bernsteinPolynomialData));
