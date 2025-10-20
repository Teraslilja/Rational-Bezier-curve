//
// (C) Matti Lehtonen 2023
//

/**
 *  @file bernstein_polynomials.hpp This file contains module interface for
 * calculation of Bernstein's polynomial factors
 */

#ifndef BERNSTEIN_POLYNOMIALS_H
#define BERNSTEIN_POLYNOMIALS_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <alloca.h>

namespace curve::bezier::utilities {

namespace factorial {
//! Maximum allowed n for factorial
static constexpr std::size_t const maximumAllowedN =
  []() consteval noexcept -> std::size_t {
  constexpr long double const availableBits =
    std::numeric_limits<std::size_t>::digits;
  long double sum = 0;
  std::size_t n = 1;
  do {
    ++n;
    sum += std::log2(static_cast<long double>(n));
  } while (sum < availableBits);
  return n - 1;
}();

/**
 *  @brief Return factorial n!
 *
 *  @param n integer value
 *  @return Factorial of n
 *  @throw std::domain:error, if n is too large
 */
[[nodiscard]] static constexpr std::size_t
factorial(std::size_t const n)
{
  if (n <= maximumAllowedN) [[likely]] {
    auto constexpr impl =
      [](std::size_t const n) constexpr noexcept -> std::size_t {
      switch (n) {
        case 0:
        case 1:
          return 1;

        case 2:
          return 2;

        [[likely]] default:
          std::size_t product = 2;
          for (std::size_t i = 3; i <= n; ++i) {
            product *= i;
          }
          return product;
      }
    };

    return impl(n);
  } else {
    throw std::domain_error("Argument n is too large!");
  }
}
}

namespace falling_factorial {
//! Maximum allowed n for falling factorial
static constexpr std::size_t const maximumAllowedN =
  []() consteval -> std::size_t {
  constexpr long double const availableBits =
    std::numeric_limits<std::size_t>::digits;

  std::size_t prevN = 2;
  do {
    std::size_t const n = prevN + 1;
    std::size_t const k = n >> 1;
    long double const bits = [=]() constexpr -> long double {
      long double sum = 0;
      for (std::size_t i = n - k + 1; i <= n; ++i) {
        sum += std::log2(static_cast<long double>(i));
      }
      return sum;
    }();
    if (bits < availableBits) {
      ++prevN;
    } else {
      return prevN;
    }
  } while (true);
}();

/**
 *  @brief Return falling factorial
 * \f$n^{\underline{k}}=n\cdot\left(n-1\right)\cdot\ldots\cdot\left(n-\left(k-1\right)\right)\f$
 *
 *  @param n integer value, n >= k >= 0
 *  @param k integer value, n >= k >= 0
 *  @return falling factorial
 */
[[nodiscard]] static constexpr std::size_t
fallingFactorial(std::size_t const n, std::size_t const k) noexcept
{
  assert(k <= n);
  if ((n == 0) || (k == 0)) {
    return 1;
  } else {
    std::size_t prod = (n - k) + 1;
    for (size_t i = prod + 1; i <= n; ++i) {
      prod *= i;
    }
    return prod;
  }
}
}

namespace binomial {
//! Implmentation interface
//! Contract: 0 <= k <= n/2, n <= maximumAllowedN
using binomialFunc = std::size_t (*)(std::size_t const,
                                     std::size_t const) noexcept;

namespace naive {
/**
 *  @brief Calculate binomial coeffient
 * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as naive
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= 2*k >= 0
 */
binomialFunc constexpr impl =
  [](std::size_t const n,
     std::size_t const k) constexpr noexcept -> std::size_t {
  std::size_t const nF = factorial::factorial(n);
  std::size_t const nkF = factorial::factorial(n - k);
  std::size_t const kF = factorial::factorial(k);

  return (nF / nkF) / kF;
};

//! Maximum allowed n for implementation
static constexpr std::size_t const maximumAllowedN = factorial::maximumAllowedN;

/**
 *  @brief Calculate binomial coeffient
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
 *  @throw std::domain_error, if k > n
 *  @throw std::domain_error, if n is too large
 */
[[nodiscard]] static constexpr std::size_t
binomial(std::size_t const n, std::size_t k)
{
  assert(k <= n);
  if (n <= maximumAllowedN) [[likely]] {
    return impl(n, std::min(k, n - k));
  } else {
    throw std::domain_error("Argument n is too large!");
  }
}
}

namespace falling_factorial {
/**
 *  @brief Calculate binomial coeffient
 * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
 *         \f$\binom{n}{k}=\begin{cases}\frac{n^{k}}{k!} &
 * \mathrm{if}\:k\leq\frac{n}{2}\\\frac{n^{\left(n-k\right)}}{\left(n-k\right)!}
 * & \mathrm{if}\:k>\frac{n}{2}\end{cases}\f$ or falling factorial
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= 2*k >= 0
 */
binomialFunc constexpr impl =
  [](std::size_t const n,
     std::size_t const k) constexpr noexcept -> std::size_t {
  std::size_t const n2k =
    curve::bezier::utilities::falling_factorial::fallingFactorial(n, k);
  std::size_t const kFact = factorial::factorial(k);

  return n2k / kFact;
};

//! Maximum allowed n for implementation
static constexpr std::size_t const maximumAllowedN =
  curve::bezier::utilities::falling_factorial::maximumAllowedN;

/**
 *  @brief Calculate binomial coeffient
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
 *  @throw std::domain_error, if n is too large
 */
[[nodiscard]] static constexpr std::size_t
binomial(std::size_t const n, std::size_t const k)
{
  assert(k <= n);
  if (n <= maximumAllowedN) [[likely]] {
    return impl(n, std::min(k, n - k));
  } else {
    throw std::domain_error("Argument n is too large!");
  }
}
}

namespace multiplication_without_recursion {
/**
 *  @brief Calculate binomial coeffient
 * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
 *         \f$\binom{n}{k}=n*\binom{n-1}{k-1}/k\f$ without recursion
 *
 * \f$\binom{6}{3}=20:\f$
 *   1
 *   1   1
 *   1   2    1
 *  _1_  3    3    1                       (3 0) = 1
 *   1  _4_   6    4   1                 4*(3 0)/1 = 4*1/1 = 4
 *   1   5  _10_  10   5 1               5*(4 1)/2 = 5*4/2 = 10
 *   1   6   15  _20_ 15 6 1 6th row (n) 6*(5 2)/3 = 6*10/3 = 20
 *                 ^ 3nd column (k)
 *
 *  @paramf n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= 2*k >= 0
 */
binomialFunc constexpr impl = [](std::size_t const n,
                                 std::size_t const k) noexcept -> std::size_t {
  switch (k) {
    case 0u:
      return 1u;

    case 1u:
      return n;

    case 2u:
      return (n * (n - 1u)) >> 1u;

    [[likely]] default:
      std::size_t diagonal = (n * (n - 1)) >> 1u;
      for (std::size_t i_k = 3, i_n = n - k; i_k < k; ++i_k, ++i_n) {
        diagonal = i_n * (diagonal / i_k);
      }

      return diagonal;
  }
};

//! Maximum allowed n for implementation
static constexpr std::size_t const maximumAllowedN =
  []() consteval noexcept -> std::size_t {
  constexpr long double const availableBits =
    std::numeric_limits<std::size_t>::digits;
  std::size_t n = 1;
  std::size_t prevBinomial = impl(n, n >> 1u);

  std::size_t bits;
  do {
    ++n;
    std::size_t const k = n >> 1u;
    std::size_t const binomial = impl(n, n >> 1u);
    bits = std::log2(static_cast<long double>(n)) +
           std::log2(static_cast<long double>(prevBinomial)) -
           std::log2(static_cast<long double>(k));
    prevBinomial = binomial;
  } while (bits < availableBits);

  return n - 1;
}();

/**
 *  @brief Calculate binomial coeffient
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
 *  @throw std::domain_error, if k > n
 *  @throw std::domain_error, if n is too large
 */
[[nodiscard]] static inline std::size_t
binomial(std::size_t const n, std::size_t k)
{
  assert(k <= n);
  if (n <= maximumAllowedN) [[likely]] {
    return impl(n, std::min(k, n - k));
  } else {
    throw std::domain_error("Argument n is too large!");
  }
}
}

namespace multiplication_with_recursion {
/**
 *  @brief Calculate binomial coeffient
 * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
 *         \f$\binom{n}{k}=\frac{n}{k}*\binom{n-1}{k-1}\f$ or recursive
 * with multiplications
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= 2*k >= 0
 */
static constexpr std::size_t
impl(std::size_t const n, std::size_t const k) noexcept
{
  switch (k) {
    case 0u:
      return 1u;

    case 1u:
      return n;

    case 2u:
      return (n * (n - 1u)) >> 1u;

    [[likely]] default:
      return n * (impl(n - 1u, k - 1u) / k);
  }
}

//! Maximum allowed n for implementation
static constexpr std::size_t const maximumAllowedN =
  []() consteval noexcept -> std::size_t {
  constexpr long double const availableBits =
    std::numeric_limits<std::size_t>::digits;
  std::size_t n = 1;
  std::size_t prevBinomial = impl(n, n >> 1u);

  std::size_t bits;
  do {
    ++n;
    std::size_t const k = n >> 1u;
    std::size_t const binomial = impl(n, n >> 1u);
    bits = std::log2(static_cast<long double>(n)) +
           std::log2(static_cast<long double>(prevBinomial)) -
           std::log2(static_cast<long double>(k));
    prevBinomial = binomial;
  } while (bits < availableBits);

  return n - 1;
}();

/**
 *  @brief Calculate binomial coeffient
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
 *  @throw std::domain_error, if n is too large
 */
[[nodiscard]] static constexpr std::size_t
binomial(std::size_t const n, std::size_t const k)
{
  assert(k <= n);
  if (n <= maximumAllowedN) [[likely]] {
    return impl(n, std::min(k, n - k));
  } else {
    throw std::domain_error("Argument n is too large!");
  }
}
}

namespace sum_without_recursion {
/**
 *  @brief Calculate binomial coeffient
 * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
 *         \f$\binom{n}{k}=\binom{n-1}{k-1}+\binom{n-1}{k}\f$ with in place
 *
 * vector \f$\binom{6}{3}=20:\f$
 *   1
 *   1  1
 *   1 _2_   1
 *   1 _3_  _3_   1
 *   1 _4_  _6_  _4_  1
 *   1  5  _10_ _10_  5 1
 *   1  6   15  _20_ 15 6 1 6th row (n)
 *                ^ 3nd column (k)
 *
 *  @paramf n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= 2*k >= 0
 */
binomialFunc const impl = [](std::size_t const n,
                             std::size_t const k) noexcept -> std::size_t {
  if ((n == 0) || (k == 0)) {
    return 1;
  } else if (k == 1) {
    return n;
  } else {
    constexpr std::size_t const maxPossibleN =
      multiplication_without_recursion::maximumAllowedN + 1;
    std::size_t* const columns = reinterpret_cast<std::size_t*>(
      alloca(sizeof(std::size_t) * (maxPossibleN << 1u)));
    std::size_t* prev = columns;
    std::size_t* curr = prev + n;
    std::size_t const nk = n - k + 1;

    // Column k=0, rows 0..n-k
    for (size_t j_n = 0; j_n < nk; ++j_n) {
      prev[j_n] = 1;
    }

    // Columns 1..k
    for (std::size_t i_k = 1; i_k <= k; ++i_k) {
      // Row k
      curr[i_k] = 1;

      // Rows k+1..(n-k) + i_k
      std::size_t const lastRow = nk + i_k;
      for (std::size_t j_n = i_k + 1; j_n < lastRow; ++j_n) {
        curr[j_n] = curr[j_n - 1] + prev[j_n - 1];
      }

      std::swap(curr, prev);
    }

    return curr[n - 1] + prev[n - 1];
  }
};

//! Maximum allowed n for implementation
static std::size_t const maximumAllowedN =
  []() constexpr noexcept -> std::size_t {
  std::size_t prevN = 2;
  std::size_t prevB = impl(prevN, prevN >> 1u);
  do {
    std::size_t const n = prevN + 1;
    std::size_t const b = impl(n, n >> 1u);
    if (b > prevB) {
      prevB = b;
      ++prevN;
    } else {
      return prevN;
    }
  } while (true);
}();

/**
 *  @brief Calculate binomial coeffient
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
 *  @throw std::domain_error, if k > n
 *  @throw std::domain_error, if n is too large
 */
[[nodiscard]] static inline std::size_t
binomial(std::size_t const n, std::size_t k)
{
  assert(k <= n);
  if (n <= maximumAllowedN) [[likely]] {
    return impl(n, std::min(k, n - k));
  } else {
    throw std::domain_error("Argument n is too large!");
  }
}
}

namespace sum_with_recursion {
/**
 *  @brief Calculate binomial coeffient
 * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
 *         \f$\binom{n}{k}=\binom{n-1}{k-1}+\binom{n-1}{k}\f$ or recursive
 * with sums
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= 2*k >= 0
 */
static constexpr std::size_t
impl(std::size_t const n, std::size_t k) noexcept
{
  switch (n) {
    case 0u:
    case 1u:
      return 1u;

    [[likely]] default:
      switch (k) {
        case 0u:
          return 1u;

        case 1u:
          return n;

        [[likely]] default:
          return impl(n - 1u, k - 1u) + impl(n - 1u, std::min(k, n - 1 - k));
      }
  }
}

//! Maximum allowed n for implementation
static std::size_t const maximumAllowedN =
  curve::bezier::utilities::binomial::sum_without_recursion::maximumAllowedN;

/**
 *  @brief Calculate binomial coeffient
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
 *  @throw std::domain_error, if n is too large
 */
[[nodiscard]] static inline std::size_t
binomial(std::size_t const n, std::size_t const k)
{
  if (n <= maximumAllowedN) [[likely]] {
    return impl(n, std::min(k, n - k));
  } else {
    throw std::domain_error("Argument n is too large!");
  }
}
}

/**
 *  @brief Calculate binomial coeffient \f$\binom{n}{k}\f$
 *
 *  @param n "k objects can be chosen from among n objects", n >= k >= 0
 *  @param k "k objects can be chosen from among n objects", n >= k >= 0
 *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
 *  @throw std::doamin_error, if k > n
 */
[[nodiscard]] static inline std::size_t
binomial(std::size_t const n, std::size_t const k)
{
  if (k <= n) [[likely]] {
    if (n < 6) {
      return multiplication_with_recursion::binomial(n, k);
    } else if (n < 9) {
      return multiplication_without_recursion::binomial(n, k);
    } else if (n <= falling_factorial::maximumAllowedN) {
      return falling_factorial::binomial(n, k);
    } else {
      return multiplication_without_recursion::binomial(n, k);
    }
  } else {
    throw std::domain_error("Argument k > n!");
  }
}
}

/**
 *  @brief Define Bernstein's polynomials and it's the first two derivatives
 * \f$B_{i,n}(u)=\binom{n}{i}u^{i}(1-u)^{n-i}\f$
 *
 *  \f$P^{th}\f$ derivate:
 *  \f$D^{p}B_{i,n}(u)=\frac{n!}{(n-p)!}\stackrel[k=\max(0,i+p-n)]{\min(i,p)}{\sum}(-1)^{k+p}\binom{p+1}{k}B_{i-k,n-p}(u)\f$
 *
 *  See e.g.
 *  Dušan Páleš, Jozef Rédl, Derivations of the Bézier Curve,
 *  Mathematics in Education, Research and Applications,
 *  htttp://dx.doi.org/10.15414/meraa.2016.02.01.01-07
 *
 *  Doha, E H; Bhrawy, A H; Saker, M A.
 *  On the Derivatives of Bernstein Polynomials: An Application for the Solution
 *  of High Even-Order Differential Equations, Boundary Value Problems; New York
 * (2011).
 */
template<typename type = float>
  requires std::is_floating_point_v<type>
struct BernsteinPolynomials
{
public:
  using real = type;

protected:
  /**
   *  @brief Return 'x' to integer power 'i'
   * \f$x^{i},x\in\mathbb{R},i\in\mathbb{N}_{0}\f$
   *
   *  @param x real value
   *  @param i ith integer power
   *  @return x to ith power
   */
  [[nodiscard]] static constexpr real toIntegerPower(real x,
                                                     std::size_t i) noexcept
  {
    real product = real(1);
    while (i > 0u) {
      if (i & 1u) {
        product *= x;
      }
      x *= x;
      i >>= 1u;
    }
    return product;
  }

public:
  /**
   *  @brief Calculate Bernstein's polynomials
   * \f$B_{i,n}(u)=\binom{n}{i}u^{i}(1-u)^{n-i}\f$
   *
   *  @param i 0 <= i <= n
   *  @param n nth degree
   *  @param u in range [0;1]
   *  @return \f$B_{i,n}(u)\f$
   *  @return zero, if i & n are not valid
   */
  [[nodiscard]] static constexpr real B(int const i,
                                        int const n,
                                        real const u) noexcept
  {
    if ((i > n) || (n < 0) || (i < 0)) {
      return real(0);
    }
    real const b = real(binomial::binomial(n, i)) * toIntegerPower(u, i) *
                   toIntegerPower(real(1) - u, n - i);
    return b;
  }

  /**
   *  @brief Calculate the first derivate of Bernstein's polynomials
   * \f$B_{i,n}(u)=\binom{n}{i}u^{i}(1-u)^{n-i}\f$
   *
   *  @param i 0 <= i <= n
   *  @param n nth degree
   *  @param u in range [0;1]
   *  @return \f$B_{i,n}(u)\f$
   *  @return zero, if i & n are not valid
   */
  [[nodiscard]] static constexpr real dB(int const i,
                                         int const n,
                                         real const u) noexcept
  {
    real const b = real(n) * (B(i - 1, n - 1, u) - B(i, n - 1, u));
    return b;
  }

  /**
   *  @brief Calculate the second derivate of Bernstein's polynomials
   * \f$B_{i,n}(u)=\binom{n}{i}u^{i}(1-u)^{n-i}\f$
   *
   *  @param i 0 <= i <= n
   *  @param n nth degree
   *  @param u in range [0;1]
   *  @return \f$B_{i,n}(u)\f$
   *  @return zero, if i & n are not valid
   */
  [[nodiscard]] static constexpr real d2B(int const i,
                                          int const n,
                                          real const u) noexcept
  {
    real const b =
      real(n) * real(n - 1) *
      (B(i, n - 2, u) - real(2) * B(i - 1, n - 2, u) + B(i - 2, n - 2, u));
    return b;
  }
};

} // namespace curve::bezier::utilities

#endif
