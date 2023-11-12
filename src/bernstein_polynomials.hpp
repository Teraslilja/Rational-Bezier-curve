//
// (C) Matti Lehtonen 2023
//

/**
 *  @file bernstein_polynomials.hpp This file contains implementation to
 * calculate Bernstein's polynomial factors
 */

#ifndef BERNSTEIN_POLYNOMIALS_H
#define BERNSTEIN_POLYNOMIALS_H

#include <cassert>
#include <type_traits>

namespace curve {
namespace bezier {
namespace internal {

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
 * of High Even-Order Differential Equations, Boundary Value Problems; New York (2011).
 */
template <typename type = float>
requires std::is_floating_point_v<type>
struct BernsteinPolynomials {
  using real = type;

protected:
  /**
   *  @brief Calculate binomial coeffient \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @return 0 otherwise
   */
  [[nodiscard]] static inline constexpr std::size_t binomial(std::size_t const n, std::size_t k) noexcept {
    if (k > n) {
      return 0u;
    } else {
      k = (k < (n >> 1u)) ? k : (n - k);
      if (k == 0u) {
        return 1u;
      } else {
        if (k == 1u) {
          return n;
        } else {
          if (k == 2u) {
            return (n * (n - 1u)) >> 1u;
          } else {
            return (binomial(n - 1u, k - 1u) * n) / k;
          }
        }
      }
    }
  }

  /**
   *  @brief Return 'x' to integer power 'i' \f$x^{i},x\in\mathbb{R},i\in\mathbb{N}_{0}\f$
   *
   *  @param x real value
   *  @param i integer power
   *  @return x to ith power
   */
  [[nodiscard]] static inline constexpr real toIntegerPower(real x, std::size_t i) noexcept {
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
   *  @brief Calculate Bernstein's polynomials \f$B_{i,n}(u)=\binom{n}{i}u^{i}(1-u)^{n-i}\f$
   *
   *  @param i 0 <= i <= n
   *  @param n nth degree
   *  @param u in range [0;1]
   *  @return \f$B_{i,n}(u)\f$
   *  @return zero, if i & n are not valid
   */
  [[nodiscard]] static inline constexpr real B(int const i, int const n, real const u) noexcept {
    if ((i > n) || (n < 0) || (i < 0)) {
      return real(0);
    }
    real const b = real(binomial(n, i)) * toIntegerPower(u, i) * toIntegerPower(real(1) - u, n - i);
    return b;
  }

  /**
   *  @brief Calculate the first derivate of Bernstein's polynomials \f$B_{i,n}(u)=\binom{n}{i}u^{i}(1-u)^{n-i}\f$
   *
   *  @param i 0 <= i <= n
   *  @param n nth degree
   *  @param u in range [0;1]
   *  @return \f$B_{i,n}(u)\f$
   *  @return zero, if i & n are not valid
   */
  [[nodiscard]] static inline constexpr real dB(int const i, int const n, real const u) noexcept {
    real const b = real(n) * (B(i - 1, n - 1, u) - B(i, n - 1, u));
    return b;
  }

  /**
   *  @brief Calculate the second derivate of Bernstein's polynomials \f$B_{i,n}(u)=\binom{n}{i}u^{i}(1-u)^{n-i}\f$
   *
   *  @param i 0 <= i <= n
   *  @param n nth degree
   *  @param u in range [0;1]
   *  @return \f$B_{i,n}(u)\f$
   *  @return zero, if i & n are not valid
   */
  [[nodiscard]] static inline constexpr real d2B(int const i, int const n, real const u) noexcept {
    real const b = real(n) * real(n - 1) * (B(i, n - 2, u) - real(2) * B(i - 1, n - 2, u) + B(i - 2, n - 2, u));
    return b;
  }
};
} // namespace internal
} // namespace bezier
} // namespace curve
#endif
