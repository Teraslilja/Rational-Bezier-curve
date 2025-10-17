//
// (C) Matti Lehtonen 2023
//

#include "bernstein_polynomials.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <format>

#include <cmath>

using real = float;

struct Testing : public curve::bezier::utilities::BernsteinPolynomials<real> {
  using Base = curve::bezier::utilities::BernsteinPolynomials<real>;

  [[nodiscard]] static constexpr std::size_t factorial(std::size_t const n) noexcept
  {
      return Base::factorial(n);
  }

  [[nodiscard]] static constexpr std::size_t fallingFactorial(std::size_t const n, std::size_t const k) noexcept
  {
      return Base::fallingFactorial(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialRecursiveMult(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialRecursiveMult(n, k);
  }

  [[nodiscard]] static inline std::size_t binomialRecursiveSum(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialRecursiveSum(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialFallingFactorial(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialFallingFactorial(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialNaive(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialNaive(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialInPlace(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialInPlace(n, k);
  }

  [[nodiscard]] static inline constexpr real toIntegerPower(real x, std::size_t i) noexcept {
    return Base::toIntegerPower(x, i);
  }
};

static constexpr std::size_t const maxFactorial_N = [] () consteval -> std::size_t
    {
        std::size_t prevN = 2;
        std::size_t prevB = Testing::factorial(prevN);
        do
        {
            std::size_t const n = prevN +1;
            std::size_t const b = Testing::factorial(n);
            if( (b/prevB) == n )
            {
                prevB = b;
                ++prevN;
            }
            else
            {
                return prevN;
            }
        }
        while(true);
    }();

static constexpr std::size_t const maxFallingFactorial_N = [] () consteval -> std::size_t
    {
        std::size_t prevN = 2;
        do
        {
            std::size_t const n = prevN +1;
            std::size_t const k = n >> 1;
            //std::size_t const b = Testing::fallingFactorial(n, k);
            long double const bits = [=]() constexpr -> long double
            {
                long double sum = 0;
                for( std::size_t i = n-k +1; i <=n; ++i )
                {
                    sum += std::log2(i);
                }
                return sum;
            }();
            if( bits < 64 )
            {
                ++prevN;
            }
            else
            {
                return prevN;
            }
        }
        while(true);
    }();

static constexpr std::size_t const maxbinomialFallingFactorial_N = maxFallingFactorial_N;

static constexpr std::size_t const maxNaive_N = maxFactorial_N;

static constexpr std::size_t const maxRecursiveMult_N = [] () consteval -> std::size_t
    {
        std::size_t prevN = 2;
        std::size_t prevB = Testing::binomialRecursiveMult(prevN, prevN>> 1u);
        do
        {
            std::size_t const n = prevN +1;
            std::size_t const b = Testing::binomialRecursiveMult(n, n >> 1u);
            //std::cerr << n << '\t' << b << std::endl;
            if( b > prevB )
            {
                prevB = b;
                ++prevN;
            }
            else
            {
                return prevN;
            }
        }
        while(true);
    }();

/*
static std::size_t const maxRecursiveSum_N = [] () -> std::size_t
    {
        std::size_t prevN = 2;
        std::size_t prevB = Testing::binomialRecursiveSum(prevN, prevN>> 1u);
        do
        {
            std::size_t const n = prevN +1;
            std::size_t const b = Testing::binomialRecursiveSum(n, n >> 1u);
            //std::cerr << n << '\t' << b << std::endl;
            if( b > prevB )
            {
                prevB = b;
                ++prevN;
            }
            else
            {
                return prevN;
            }
        }
        while(true);
    }();
*/

static std::size_t const maxInPlace_N = [] () -> std::size_t
    {
        std::size_t prevN = 2;
        std::size_t prevB = Testing::binomialInPlace(prevN, prevN>> 1u);
        do
        {
            std::size_t const n = prevN +1;
            std::size_t const b = Testing::binomialInPlace(n, n >> 1u);
            //std::cerr << n << '\t' << b << std::endl;
            if( b > prevB )
            {
                prevB = b;
                ++prevN;
            }
            else
            {
                return prevN;
            }
        }
        while(true);
    }();

struct FactorialData {
  std::size_t const n;

  std::size_t const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, FactorialData const &data) {
    out << "{";
    out << data.n << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class FactorialTests : public ::testing::TestWithParam<FactorialData> {};

TEST_P(FactorialTests, factorial_correctResult) {
  auto const &params = GetParam();

  real const result = Testing::factorial(params.n);

  EXPECT_EQ(result, params.expectedResult);
}

static constexpr FactorialData factorialData[]{
    {0, 1},
    {1, 1},
    {2, 1*2},
    {3, 1*2*3},
    {4, 1*2*3*4},
    {5, 1*2*3*4*5},
    {6, 1*2*3*4*5*6},
    {7, 1*2*3*4*5*6*7},
    {20, 1ULL*2*3*4*5*6*7*8*9*10*11*12*13*14*15*16*17*18*19*20}, //< The largest factorial for 64 bits
};

INSTANTIATE_TEST_SUITE_P(Fixture, FactorialTests, ::testing::ValuesIn(factorialData));

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

TEST(BinomialTests, maximum_N)
{
    [[maybe_unused]] std::size_t maxVal;

    std::cerr << "maxFactorial_N = " << maxFactorial_N << std::endl;
    maxVal = Testing::factorial(maxFactorial_N);
    //std::cerr << "maxVal = " << maxVal << std::endl;

    std::cerr << "maxFallingFactorial_N = " << maxFallingFactorial_N << std::endl;
    maxVal = Testing::fallingFactorial(maxFallingFactorial_N, maxFallingFactorial_N>>1u);
    //std::cerr << "maxVal = " << maxVal << std::endl<<std::endl;
    
    std::cerr << "maxbinomialFallingFactorial_N = " << maxbinomialFallingFactorial_N << std::endl;
    maxVal = Testing::binomialFallingFactorial(maxbinomialFallingFactorial_N, maxbinomialFallingFactorial_N >> 1u);
    //std::cerr << "maxVal = " << maxVal << std::endl;

    std::cerr << "maxNaive_N = " << maxNaive_N << std::endl;
    maxVal = Testing::binomialNaive(maxNaive_N, maxNaive_N >> 1u);
    //std::cerr << "maxVal = " << maxVal << std::endl;

    std::cerr << "maxRecursiveMult_N = " << maxRecursiveMult_N << std::endl;
    maxVal = Testing::binomialRecursiveMult(maxRecursiveMult_N, maxRecursiveMult_N >> 1u);
    //std::cerr << "maxVal = " << maxVal << std::endl;
    
    //std::cerr << "maxRecursiveSum_N = " << maxRecursiveSum_N << std::endl;
    //maxVal = Testing::binomialRecursiveSum(maxRecursiveSum_N, maxRecursiveSum_N >> 1u);
    //std::cerr << "maxVal = " << maxVal << std::endl;
    
    std::cerr << "maxInPlace_N = " << maxInPlace_N << std::endl;
    maxVal = Testing::binomialInPlace(maxInPlace_N, maxInPlace_N >> 1u);
    //std::cerr << "maxVal = " << maxVal << std::endl;
}

class BinomialTests : public ::testing::TestWithParam<BinomialData> {};

TEST_P(BinomialTests, binomialRecursiveMult_correctResult) {
  auto const &params = GetParam();

  std::size_t const result = Testing::binomialRecursiveMult(params.n, params.k);

  EXPECT_EQ(result, params.expectedResult);
}

TEST_P(BinomialTests, binomialRecursiveSum_correctResult) {
  auto const &params = GetParam();

  std::size_t const result = Testing::binomialRecursiveSum(params.n, params.k);

  EXPECT_EQ(result, params.expectedResult);
}

TEST_P(BinomialTests, binomialFallingFactorial_correctResult) {
  auto const &params = GetParam();

  std::size_t const result = Testing::binomialFallingFactorial(params.n, params.k);

  EXPECT_EQ(result, params.expectedResult);
}

TEST_P(BinomialTests, binomialNaive_correctResult) {
  auto const &params = GetParam();

  std::size_t const result = Testing::binomialNaive(params.n, params.k);

  EXPECT_EQ(result, params.expectedResult);
}

TEST_P(BinomialTests, binomialInPlace_correctResult) {
  auto const &params = GetParam();

  std::size_t const result = Testing::binomialInPlace(params.n, params.k);

  EXPECT_EQ(result, params.expectedResult);
}

//      1
//     1  1
//    1  2  1
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
    {6u, 7u, 0u},
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


struct FallingFactorialData {
  std::size_t const n;
  std::size_t const k;

  std::size_t const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, FallingFactorialData const &data) {
    out << "{";
    out << data.n << "," << data.k << "," << data.expectedResult;
    out << "}";
    return out;
  }
};

class FallingFactorialTests: public ::testing::TestWithParam<FallingFactorialData> {};

TEST_P(FallingFactorialTests, fallingFactorial_correctResult) {
  auto const &params = GetParam();

    std::size_t const result = Testing::fallingFactorial( params.n, params.k );

    EXPECT_EQ(result, params.expectedResult);
}

static constexpr FallingFactorialData fallingFactorialData[]{
    { 0, 0, 1},
    { 1, 0, 1}, { 1, 1, 1},
    { 2, 0, 1}, { 2, 1, 2}, { 2, 2, 2*1},
    { 3, 0, 1}, { 3, 1, 3}, { 3, 2, 3*2}, { 3, 3, 3*2*1},
    { 4, 0, 1}, { 4, 1, 4}, { 4, 2, 4*3}, { 4, 3, 4*3*2}, { 4, 4, 4*3*2*1},
    { 5, 0, 1}, { 5, 1, 5}, { 5, 2, 5*4}, { 5, 3, 5*4*3}, { 5, 4, 5*4*3*2}, { 5, 5, 5*4*3*2*1},
    { 6, 0, 1}, { 6, 1, 6}, { 6, 2, 6*5}, { 6, 3, 6*5*4}, { 6, 4, 6*5*4*3}, { 6, 5, 6*5*4*3*2}, { 6, 6, 6*5*4*3*2*1},
    { 28, 14, 28ULL*27*26*25*24*23*22*21*20*19*18*17*16*15},
};

INSTANTIATE_TEST_SUITE_P(Fixture, FallingFactorialTests, ::testing::ValuesIn(fallingFactorialData));


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
