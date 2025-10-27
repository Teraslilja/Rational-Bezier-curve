//
// (C) Matti Lehtonen 2023
//

/**
 *  @file bernstein_polynomials.hpp This file contains module interface for
 * calculation of Bernstein's polynomial factors
 */

#ifndef BERNSTEIN_POLYNOMIALS_H
#define BERNSTEIN_POLYNOMIALS_H

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <alloca.h>

namespace curve::bezier::utilities {

template<typename UINT = std::size_t, typename INDEX = std::size_t>
  requires std::is_integral_v<UINT> && std::is_integral_v<INDEX> &&
           std::is_unsigned_v<UINT> && std::is_unsigned_v<INDEX>
struct Factorial
{
public:
  using Index = INDEX;
  using Integer = UINT;

private:
  static constexpr Integer impl(Index const n) noexcept
  {
    switch (n) {
      case 0:
      case 1:
        return Integer(1);

      case 2:
        return Integer(2);

      [[likely]] default:
        Integer product = Integer(2);
        for (Index i = 3; i <= n; ++i) {
          product *= Integer(i);
        }
        return product;
    }
  }

protected:
  //! Maximum allowed n for factorial
  static constexpr Index const maximumAllowedN =
    []() consteval noexcept -> Index {
    constexpr long double const availableBits =
      std::numeric_limits<Integer>::digits;
    long double sum = 0;
    Index n = 1;
    do {
      ++n;
      sum += std::log2(static_cast<long double>(n));
    } while (sum < availableBits);
    return n - 1;
  }();

public:
  /**
   *  @brief Return factorial n!
   *
   *  @param n integer value
   *  @return Factorial of n
   *  @throw std::domain:error, if n is too large
   */
  [[nodiscard]] static constexpr Integer factorial(Index const n)
  {
    if (n > maximumAllowedN) [[unlikely]] {
      throw std::domain_error("Argument n is too large!");
    }
    return impl(n);
  }
};

template<typename UINT = std::size_t, typename INDEX = std::size_t>
  requires std::is_integral_v<UINT> && std::is_integral_v<INDEX> &&
           std::is_unsigned_v<UINT> && std::is_unsigned_v<INDEX>
struct FallingFactorial
{
public:
  using Index = INDEX;
  using Integer = UINT;

private:
  static constexpr Integer impl(Index const n, Index const k) noexcept
  {
    if ((n == Index(0)) || (k == Index(0))) {
      return Integer(1);
    } else {
      Integer prod = (n - k) + 1;
      for (Index i = Index(prod) + 1; i <= n; ++i) {
        prod *= Index(i);
      }
      return prod;
    }
  }

protected:
  //! Maximum allowed n for falling factorial
  static constexpr Index const maximumAllowedN = []() consteval -> Index {
    constexpr long double const availableBits =
      std::numeric_limits<Integer>::digits;

    Index prevN = 2;
    do {
      Index const n = prevN + 1;
      Index const k = n >> 1;
      long double const bits = [=]() constexpr -> long double {
        long double sum = 0;
        for (Index i = n - k + 1; i <= n; ++i) {
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

public:
  /**
   *  @brief Return falling factorial
   * \f$n^{\underline{k}}=n\cdot\left(n-1\right)\cdot\ldots\cdot\left(n-\left(k-1\right)\right)\f$
   *
   *  @param n integer value, n >= k >= 0
   *  @param k integer value, n >= k >= 0
   *  @return falling factorial
   *  @throw std::domain_error, if k > n
   *  @throw std::domain:error, if n is too large
   */
  [[nodiscard]] static constexpr Integer fallingFactorial(Index const n,
                                                          Index const k)
  {
    if (n > maximumAllowedN) [[unlikely]] {
      throw std::domain_error("Argument n is too large!");
    }
    if (k > n) [[unlikely]] {
      throw std::domain_error("Argument k > argument n!");
    }
    return impl(n, k);
  }
};

namespace binomial {
struct Naive : private Factorial<std::size_t, std::size_t>
{
private:
  using Base = Factorial<std::size_t, std::size_t>;

  /**
   *  @brief Calculate binomial coeffient
   * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as naive
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= 2*k >= 0
   */
  static constexpr std::size_t impl(std::size_t const n,
                                    std::size_t const k) noexcept
  {
    std::size_t const nF = factorial(n);
    std::size_t const nkF = factorial(n - k);
    std::size_t const kF = factorial(k);

    return (nF / nkF) / kF;
  }

public:
  //! Maximum allowed n for implementation
  using Base::maximumAllowedN;

  /**
   *  @brief Calculate binomial coeffient
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @throw std::domain_error, if k > n
   *  @throw std::domain_error, if n is too large
   */
  [[nodiscard]] static constexpr std::size_t binomial(std::size_t const n,
                                                      std::size_t k)
  {
    if (n > maximumAllowedN) [[unlikely]] {
      throw std::domain_error("Argument n is too large!");
    }
    return impl(n, std::min(k, n - k));
  }
};

struct FallingFactorial
  : private Factorial<std::size_t, std::size_t>
  , private curve::bezier::utilities::FallingFactorial<std::size_t, std::size_t>
{
private:
  using BaseF = Factorial<std::size_t, std::size_t>;
  using BaseFF =
    curve::bezier::utilities::FallingFactorial<std::size_t, std::size_t>;

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
  static constexpr std::size_t impl(std::size_t const n,
                                    std::size_t const k) noexcept
  {
    std::size_t const n2k = fallingFactorial(n, k);
    std::size_t const kFact = factorial(k);

    return n2k / kFact;
  }

public:
  //! Maximum allowed n for implementation
  static constexpr std::size_t const maximumAllowedN = BaseFF::maximumAllowedN;

  //! Maximum allowed k for implementation
  static constexpr std::size_t const maximumAllowedK = BaseF::maximumAllowedN;

  /**
   *  @brief Calculate binomial coeffient
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @throw std::domain_error, if n is too large
   *  @throw std::domain_error, if k is too large
   */
  [[nodiscard]] static constexpr std::size_t binomial(std::size_t const n,
                                                      std::size_t const k)
  {
    if (n > maximumAllowedN) [[unlikely]] {
      throw std::domain_error("Argument n is too large!");
    }
    if (k > maximumAllowedK) [[unlikely]] {
      throw std::domain_error("Argument k is too large!");
    }

    return impl(n, std::min(k, n - k));
  }
};

struct MultiplicationWithoutRecursion
{
private:
  /**
   *  @brief Calculate binomial coeffient
   * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
   *         \f$\binom{n}{k}=n*\binom{n-1}{k-1}/k\f$ without recursion
   *
   *
   * \f$\binom{4}{0} =  1 =        (1)/(0!)\f$
   * \f$\binom{5}{1} =  5 =        (5)/(1)       = 5*\binom{4}{0}/1\f$
   * \f$\binom{6}{2} = 15 =      (6*5)/(2*1)     = 6*\binom{5}{1}/2\f$
   * \f$\binom{7}{3} = 35 =    (7*6*5)/(3*2*1)   = 7*\binom{6}{2}/3\f$
   * \f$\binom{8}{4} = 70 =  (8*7*6*5)/(4*3*2*1) = 8*\binom{7}{3}/4\f$
   *
   *  @paramf n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= 2*k >= 0
   */
  static constexpr std::size_t impl(std::size_t const n,
                                    std::size_t const k) noexcept
  {
    switch (k) {
      case 0u:
        return 1u;

      case 1u:
        return n;

      case 2u:
        return (n * (n - 1u)) >> 1u;

      [[likely]] default:
        std::size_t i_k = 2;
        std::size_t i_n = n - (k - i_k);
        std::size_t diagonal = (i_n * (i_n - 1u)) >> 1u;
        for (++i_k, ++i_n; i_k <= k; ++i_k, ++i_n) {
          diagonal = (i_n * diagonal) / i_k;
        }

        return diagonal;
    }
  };

public:
  //! Maximum allowed n for implementation
  static constexpr std::size_t const maximumAllowedN =
    []() consteval noexcept -> std::size_t {
    auto constexpr I =
      [](std::size_t const n,
         std::size_t const k) consteval noexcept -> std::size_t {
      std::size_t i_k = 2;
      std::size_t i_n = n - (k - i_k);
      std::size_t diagonal = (i_n * (i_n - 1u)) >> 1u;
      for (++i_k, ++i_n; i_k <= k; ++i_k, ++i_n) {
        diagonal = (i_n * diagonal) / i_k;
      }

      return diagonal;
    };
    std::size_t n = FallingFactorial::maximumAllowedN;
    std::size_t k = n >> 1u;
    std::size_t prev = I(n, k);
    std::size_t next;
    for (++n, k = n >> 1u; prev < (next = I(n, k)); ++n, k = n >> 1u) {
      prev = next;
    }
    return n - 1u;
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
  [[nodiscard]] static inline std::size_t binomial(std::size_t const n,
                                                   std::size_t k)
  {
    if (n > maximumAllowedN) [[unlikely]] {
      throw std::domain_error("Argument n is too large!");
    }

    return impl(n, std::min(k, n - k));
  }
};

struct MultiplicationWithRecursion
{
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
  static constexpr std::size_t impl(std::size_t const n,
                                    std::size_t const k) noexcept
  {
    switch (k) {
      case 0u:
        return 1u;

      case 1u:
        return n;

      case 2u:
        return (n * (n - 1u)) >> 1u;

      [[likely]] default:
        return (n * impl(n - 1u, k - 1u)) / k;
    }
  }

  //! Maximum allowed n for implementation
  static constexpr std::size_t const maximumAllowedN =
    MultiplicationWithoutRecursion::maximumAllowedN;

  /**
   *  @brief Calculate binomial coeffient
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @throw std::domain_error, if n is too large
   */
  [[nodiscard]] static constexpr std::size_t binomial(std::size_t const n,
                                                      std::size_t const k)
  {
    if (n > maximumAllowedN) [[unlikely]] {
      throw std::domain_error("Argument n is too large!");
    }
    return impl(n, std::min(k, n - k));
  }
};

struct SumWithoutRecursion
{
  //! Maximum allowed n for implementation
  static std::size_t const maximumAllowedN =
    []() consteval noexcept -> std::size_t {
    constexpr long double const maxBits =
      std::numeric_limits<std::size_t>::digits;
    for (std::size_t n = FallingFactorial::maximumAllowedN;; ++n) {
      std::size_t k = n >> 1u;
      std::size_t nk = n - k;

      long double bits = 0;
      for (std::size_t i_n = n; i_n > nk; --i_n) {
        bits += std::log2(static_cast<long double>(i_n));
      }
      for (std::size_t i_k = 2u; i_k <= k; ++i_k) {
        bits -= std::log2(static_cast<long double>(i_k));
      }
      if (bits >= maxBits) {
        return n - 1u;
      }
    }
  }();

  /**
   *  @brief Calculate binomial coeffient
   * \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
   *         \f$\binom{n}{k}=\binom{n-1}{k-1}+\binom{n-1}{k}\f$ with in place
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
  static constexpr std::size_t impl(std::size_t const n,
                                    std::size_t const k) noexcept
  {
    if ((n == 0u) || (k == 0u)) {
      return 1u;
    } else if (k == 1u) {
      return n;
    } else {
      std::array<std::size_t, maximumAllowedN + 1u> column1;
      std::array<std::size_t, maximumAllowedN + 1u> column2;
      std::size_t* prev = column1.data();
      std::size_t* curr = column2.data();
      std::size_t const nk = n - k + 1u;

      // Column k=0, rows 0..n-k
      for (size_t j_n = 0u; j_n < nk; ++j_n) {
        prev[j_n] = 1u;
      }

      // Columns 1..k
      for (std::size_t i_k = 1u; i_k <= k; ++i_k) {
        // Row k
        curr[i_k] = 1u;

        // Rows k+1..(n-k) + i_k
        std::size_t const lastRow = nk + i_k;
        for (std::size_t j_n = i_k + 1u; j_n < lastRow; ++j_n) {
          curr[j_n] = curr[j_n - 1u] + prev[j_n - 1u];
        }

        std::swap(curr, prev);
      }

      return curr[n - 1u] + prev[n - 1u];
    }
  }

  /**
   *  @brief Calculate binomial coeffient
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @throw std::domain_error, if n is too large
   */
  [[nodiscard]] static constexpr std::size_t binomial(std::size_t const n,
                                                      std::size_t const k)
  {
    if (n > maximumAllowedN) [[unlikely]] {
      throw std::domain_error("Argument n is too large!");
    }
    return impl(n, std::min(k, n - k));
  }
};

struct SumWithRecursion
{
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
  static constexpr std::size_t impl(std::size_t const n, std::size_t k) noexcept
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
            return impl(n - 1u, k - 1u) + impl(n - 1u, std::min(k, n - k - 1u));
        }
    }
  }

  //! Maximum allowed n for implementation, very slow ..
  static std::size_t constexpr maximumAllowedN = std::min(
    std::size_t(25),
    curve::bezier::utilities::binomial::SumWithoutRecursion::maximumAllowedN);

  /**
   *  @brief Calculate binomial coeffient
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @throw std::domain_error, if n is too large
   */
  [[nodiscard]] static constexpr std::size_t binomial(std::size_t const n,
                                                      std::size_t const k)
  {
    if (n > maximumAllowedN) [[unlikely]] {
      throw std::domain_error("Argument n is too large!");
    }
    return impl(n, std::min(k, n - k));
  }
};

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
  if (k > n) [[unlikely]] {
    throw std::domain_error("Argument k > n!");
  }
  return SumWithoutRecursion::binomial(n, k);
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
