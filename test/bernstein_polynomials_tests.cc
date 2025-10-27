//
// (C) Matti Lehtonen 2023
//

#include "bernstein_polynomials.hpp"

#include <gtest/gtest.h>

#include <format>
#include <iostream>

#include <cmath>

using real = float;

using IndexType = std::size_t;
using IntegerType = std::size_t;

using namespace curve::bezier::utilities;

struct FactorialTesting : Factorial<IntegerType, IndexType>
{
  using Base = Factorial<IntegerType, IndexType>;
  using Base::maximumAllowedN;

  [[nodiscard]] static constexpr IntegerType factorial(
    IndexType const n) noexcept
  {
    return Base::factorial(n);
  }
};

TEST(FactorialTests, maximum_N)
{
  std::cerr << "FactorialTesting::maximumAllowedN = "
            << FactorialTesting::maximumAllowedN << std::endl;
}

struct FactorialData
{
  IndexType const n;

  IntegerType const expectedResult;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         FactorialData const& data)
  {
    out << "{";
    out << data.n << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class FactorialTests : public ::testing::TestWithParam<FactorialData>
{};

TEST_P(FactorialTests, factorial_correctResult)
{
  auto const& params = GetParam();

  if (params.n <= FactorialTesting::maximumAllowedN) {
    real const result = FactorialTesting::factorial(params.n);

    EXPECT_EQ(result, params.expectedResult);
  }
}

static constexpr FactorialData factorialData[]{
  { 0, 1 },
  { 1, 1 },
  { 2, 1 * 2 },
  { 3, 1 * 2 * 3 },
  { 4, 1 * 2 * 3 * 4 },
  { 5, 1 * 2 * 3 * 4 * 5 },
  { 6, 1 * 2 * 3 * 4 * 5 * 6 },
  { 7, 1 * 2 * 3 * 4 * 5 * 6 * 7 },
  { 20,
    std::size_t(1) * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11 * 12 * 13 * 14 *
      15 * 16 * 17 * 18 * 19 * 20 }, //< The largest factorial for 64 bits
};

INSTANTIATE_TEST_SUITE_P(Fixture,
                         FactorialTests,
                         ::testing::ValuesIn(factorialData));

struct FallingFactorialTesting : FallingFactorial<IntegerType, IndexType>
{
  using Base = FallingFactorial<IntegerType, IndexType>;
  using Base::maximumAllowedN;

  [[nodiscard]] static constexpr IntegerType fallingFactorial(
    IndexType const n,
    IndexType const k) noexcept
  {
    return Base::fallingFactorial(n, k);
  }
};

struct FallingFactorialData
{
  IndexType const n;
  IndexType const k;

  IntegerType const expectedResult;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         FallingFactorialData const& data)
  {
    out << "{";
    out << data.n << "," << data.k << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

TEST(FallingFactorialTests, maximum_N)
{
  std::cerr << "FallingFactorial::maximumAllowedN = "
            << FallingFactorialTesting::maximumAllowedN << std::endl;
}

class FallingFactorialTests
  : public ::testing::TestWithParam<FallingFactorialData>
{};

TEST_P(FallingFactorialTests, fallingFactorial_correctResult)
{
  auto const& params = GetParam();

  if (params.n <= FallingFactorialTesting::maximumAllowedN) {
    real const result =
      FallingFactorialTesting::fallingFactorial(params.n, params.k);

    EXPECT_EQ(result, params.expectedResult);
  }
}

/* clang-format off */
static constexpr FallingFactorialData fallingFactorialData[]{
  { 0, 0, 1 },
  { 1, 0, 1 }, { 1, 1, 1 },
  { 2, 0, 1 }, { 2, 1, 2 }, { 2, 2, 2 * 1 },
  { 3, 0, 1 }, { 3, 1, 3 }, { 3, 2, 3 * 2 }, { 3, 3, 3 * 2 * 1 },
  { 4, 0, 1 }, { 4, 1, 4 }, { 4, 2, 4 * 3 }, { 4, 3, 4 * 3 * 2 }, { 4, 4, 4 * 3 * 2 * 1 },
  { 5, 0, 1 }, { 5, 1, 5 }, { 5, 2, 5 * 4 }, { 5, 3, 5 * 4 * 3 }, { 5, 4, 5 * 4 * 3 * 2 }, { 5, 5, 5 * 4 * 3 * 2 * 1 },
  { 6, 0, 1 }, { 6, 1, 6 }, { 6, 2, 6 * 5 }, { 6, 3, 6 * 5 * 4 }, { 6, 4, 6 * 5 * 4 * 3 }, { 6, 5, 6 * 5 * 4 * 3 * 2 }, { 6, 6, 6 * 5 * 4 * 3 * 2 * 1 },
  { 29, 14, std::size_t(29) * 28 * 27 * 26 * 25 * 24 * 23 * 22 * 21 * 20 * 19 * 18 * 17 * 16 },
};
/* clang-format on */

INSTANTIATE_TEST_SUITE_P(Fixture,
                         FallingFactorialTests,
                         ::testing::ValuesIn(fallingFactorialData));

struct BinomialTesting
{
  [[nodiscard]] static constexpr std::size_t binomialNaive(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return binomial::Naive::binomial(n, k);
  }

  [[nodiscard]] static constexpr std::size_t binomialFallingFactorial(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return binomial::FallingFactorial::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t
  binomialMultiplicationWithoutRecursion(std::size_t const n,
                                         std::size_t const k) noexcept
  {
    return binomial::MultiplicationWithoutRecursion::binomial(n, k);
  }

  [[nodiscard]] static constexpr std::size_t
  binomialMultiplicationWithRecursion(std::size_t const n,
                                      std::size_t const k) noexcept
  {
    return binomial::MultiplicationWithRecursion::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t binomialSumWithoutRecursion(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return binomial::SumWithoutRecursion::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t binomialSumWithRecursion(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return binomial::SumWithRecursion::binomial(n, k);
  }
};

struct BinomialData
{
  std::size_t const n;
  std::size_t const k;

  std::size_t const expectedResult;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         BinomialData const& data)
  {
    out << "{";
    out << data.n << "," << data.k << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

TEST(BinomialTests, maximum_N)
{
  std::cerr << "binomial::Naive::maximumAllowedN = "
            << binomial::Naive::maximumAllowedN << std::endl;

  std::cerr << "binomial::FallingFactorial::maximumAllowedN = "
            << binomial::FallingFactorial::maximumAllowedN << std::endl;

  std::cerr << "binomial::MultiplicationWithoutRecursion::maximumAllowedN = "
            << binomial::MultiplicationWithoutRecursion::maximumAllowedN
            << std::endl;

  std::cerr << "binomial::MultiplicationWithRecursion::maximumAllowedN = "
            << binomial::MultiplicationWithRecursion::maximumAllowedN
            << std::endl;

  std::cerr << "binomial::SumWithoutRecursion::maximumAllowedN = "
            << binomial::SumWithoutRecursion::maximumAllowedN << std::endl;

  std::cerr << "binomial::SumWithRecursion::maximumAllowedN = "
            << binomial::SumWithRecursion::maximumAllowedN << std::endl;
}

class BinomialTests : public ::testing::TestWithParam<BinomialData>
{};

TEST_P(BinomialTests, binomialNaive_correctResult)
{
  auto const& params = GetParam();

  if (params.n <= binomial::Naive::maximumAllowedN) {
    std::size_t const result =
      BinomialTesting::binomialNaive(params.n, params.k);

    EXPECT_EQ(result, params.expectedResult);
  }
}

TEST_P(BinomialTests, binomialFallingFactorial_correctResult)
{
  auto const& params = GetParam();

  if (params.n <= binomial::FallingFactorial::maximumAllowedN) {
    std::size_t const result =
      BinomialTesting::binomialFallingFactorial(params.n, params.k);

    EXPECT_EQ(result, params.expectedResult);
  }
}

TEST_P(BinomialTests, binomialMultiplicationWithoutRecursion_correctResult)
{
  auto const& params = GetParam();

  if (params.n <= binomial::MultiplicationWithoutRecursion::maximumAllowedN) {
    std::size_t const result =
      BinomialTesting::binomialMultiplicationWithoutRecursion(params.n,
                                                              params.k);

    EXPECT_EQ(result, params.expectedResult);
  }
}

TEST_P(BinomialTests, binomialRecursiveMult_correctResult)
{
  auto const& params = GetParam();

  if (params.n <= binomial::MultiplicationWithRecursion::maximumAllowedN) {
    std::size_t const result =
      BinomialTesting::binomialMultiplicationWithRecursion(params.n, params.k);

    EXPECT_EQ(result, params.expectedResult);
  }
}

TEST_P(BinomialTests, binomialSumWithoutRecursion_correctResult)
{
  auto const& params = GetParam();

  if (params.n <= binomial::SumWithoutRecursion::maximumAllowedN) {
    std::size_t const result =
      BinomialTesting::binomialSumWithoutRecursion(params.n, params.k);

    EXPECT_EQ(result, params.expectedResult);
  }
}

TEST_P(BinomialTests, binomialSumWithRecursion_correctResult)
{
  auto const& params = GetParam();

  if (params.n <= binomial::SumWithRecursion::maximumAllowedN) {
    std::size_t const result =
      BinomialTesting::binomialSumWithRecursion(params.n, params.k);

    EXPECT_EQ(result, params.expectedResult);
  }
}

TEST_P(BinomialTests, binomial_correctResult)
{
  auto const& params = GetParam();

  std::size_t const result = binomial::binomial(params.n, params.k);

  EXPECT_EQ(result, params.expectedResult);
}

//      1
//     1  1
//    1  2  1
//   1 3   3 1
//  1 4  6  4 1
// 1 5 10 10 5 1
/* clang-format off */
static constexpr BinomialData binomialData[]{
  { 0u, 0u, 1u }, //
  { 1u, 0u, 1u }, { 1u, 1u, 1u }, //
  { 2u, 0u, 1u }, { 2u, 1u, 2u }, { 2u, 2u, 1u }, //
  { 3u, 0u, 1u }, { 3u, 1u, 3u }, { 3u, 2u, 3u }, { 3u, 3u, 1u }, //
  { 4u, 0u, 1u }, { 4u, 1u, 4u }, { 4u, 2u, 6u }, { 4u, 3u, 4u }, { 4u, 4u, 1u }, //
  { 5u, 0u, 1u }, { 5u, 1u, 5u }, { 5u, 2u, 10u }, { 5u, 3u, 10u }, { 5u, 4u, 5u }, { 5u, 5u, 1u }, //
  { 20u, 10u, 184756u },
  { 29u, 14u, 77558760u },
  { 62u, 31u, 465428353255261088u },
  { 67u, 33u, 14226520737620288370u },
  // { 68u, 34u, 28453041475240576740u },
  // { 69u, 34u, 56093138908331422716u },
  //{ 70u, 35u, 112186277816662845432u },
};
/* clang-format on */

INSTANTIATE_TEST_SUITE_P(Fixture,
                         BinomialTests,
                         ::testing::ValuesIn(binomialData));

struct Testing : public curve::bezier::utilities::BernsteinPolynomials<real>
{
  using Base = curve::bezier::utilities::BernsteinPolynomials<real>;

  static constexpr real toIntegerPower(real const x, std::size_t const i)
  {
    return Base::toIntegerPower(x, i);
  }
};

struct IntegerPowerData
{
  real const x;
  std::size_t const i;

  real const expectedResult;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         IntegerPowerData const& data)
  {
    out << "{";
    out << data.x << "," << data.i << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class IntegerPowerTests : public ::testing::TestWithParam<IntegerPowerData>
{};

TEST_P(IntegerPowerTests, integerPower_correctResult)
{
  auto const& params = GetParam();

  real const result = Testing::toIntegerPower(params.x, params.i);

  EXPECT_EQ(result, params.expectedResult);
}

/* clang-format off */
static constexpr IntegerPowerData integerPowerData[]{
  { real(0), 0u, real(1) }, { real(0.5), 0u, real(1) }, { real(1), 0u, real(1) },

  { real(0), 1u, real(0) }, { real(0.5), 1u, real(0.5) }, { real(1), 1u, real(1) },

  { real(0), 2u, real(0) }, { real(0.5), 2u, real(0.25) }, { real(1), 2u, real(1) },

  { real(0), 3u, real(0) }, { real(0.5), 3u, real(0.125) }, { real(1), 3u, real(1) },
};
/* clang-format on */

INSTANTIATE_TEST_SUITE_P(Fixture,
                         IntegerPowerTests,
                         ::testing::ValuesIn(integerPowerData));

struct BernsteinPolynomialData
{
  std::size_t const i;
  std::size_t const n;
  real const u;

  real const expectedResultB;
  real const expectedResultdB;

  inline friend std::ostream& operator<<(std::ostream& out,
                                         BernsteinPolynomialData const& data)
  {
    out << "{";
    out << data.i << "," << data.n << "," << data.u << ","
        << data.expectedResultB << "," << data.expectedResultdB;
    out << "}";
    return out;
  }
};

class BernsteinPolynomialTests
  : public ::testing::TestWithParam<BernsteinPolynomialData>
{};

TEST_P(BernsteinPolynomialTests, B_correctResult)
{
  auto const& params = GetParam();

  real const result = Testing::B(params.i, params.n, params.u);

  EXPECT_EQ(result, params.expectedResultB);
}

TEST_P(BernsteinPolynomialTests, dB_correctResult)
{
  auto const& params = GetParam();

  real const result = Testing::dB(params.i, params.n, params.u);

  EXPECT_EQ(result, params.expectedResultdB);
}

static constexpr BernsteinPolynomialData bernsteinPolynomialData[]{
  { 0u, 1u, real(0), real(1 * 1.0 * 1.0), real(1 * (0.0 - 1.0)) },
  { 0u, 1u, real(0.5), real(1 * 1.0 * 0.5), real(1 * (0.0 - 1.0)) },
  { 0u, 1u, real(1), real(1 * 1.0 * 0.0), real(1 * (0.0 - 1.0)) },

  { 0u, 2u, real(0), real(1 * 1.0 * 1.0), real(2 * (0.0 - 1.0)) },
  { 0u, 2u, real(0.5), real(1 * 1.0 * 0.25), real(2 * (0.0 - 0.5)) },
  { 0u, 2u, real(1), real(1 * 1.0 * 0.0), real(2 * (0.0 - 0.0)) },

  { 1u, 2u, real(0), real(2 * 0.0 * 1.0), real(2 * (1.0 - 0.0)) },
  { 1u, 2u, real(0.5), real(2 * 0.5 * 0.5), real(2 * (0.5 - 0.5)) },
  { 1u, 2u, real(1), real(2 * 1.0 * 0.0), real(2 * (0.0 - 1.0)) },

  { 1u, 3u, real(0), real(3 * 0.0 * 1.0), real(3 * (1.0 - 0.0)) },
  { 1u, 3u, real(0.5), real(3 * 0.5 * 0.25), real(3 * (0.25 - 0.5)) },
  { 1u, 3u, real(1), real(3 * 1.0 * 0.0), real(3 * (0.0 - 0.0)) },

  { 2u, 3u, real(0), real(3 * 0.0 * 1.0), real(3 * (0.0 - 0.0)) },
  { 2u, 3u, real(0.5), real(3 * 0.25 * 0.5), real(3 * (0.5 - 0.25)) },
  { 2u, 3u, real(1), real(3 * 1.0 * 0.0), real(3 * (0.0 - 1.0)) },
};

INSTANTIATE_TEST_SUITE_P(Fixture,
                         BernsteinPolynomialTests,
                         ::testing::ValuesIn(bernsteinPolynomialData));
