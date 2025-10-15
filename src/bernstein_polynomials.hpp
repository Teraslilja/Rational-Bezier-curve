//
// (C) Matti Lehtonen 2023
//

/**
 *  @file bernstein_polynomials.hpp This file contains module interface for calculation of Bernstein's polynomial
 * factors
 */

#ifndef BERNSTEIN_POLYNOMIALS_H
#define BERNSTEIN_POLYNOMIALS_H

#include <type_traits>
#include <vector>

namespace curve::bezier::utilities {

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
public:
  using real = type;

protected:
  /**
   *  @brief Calculate binomial coeffient \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
   *         \f$\binom{n}{k}=\binom{n-1}{k-1}+\binom{n-1}{k}\f$ or recursive
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @return 0 otherwise
   */
  [[nodiscard]] static constexpr std::size_t binomialRecursive(std::size_t const n, std::size_t const k) noexcept {
    auto constexpr bin = [](std::size_t const n, std::size_t const k) -> std::size_t {
      auto constexpr impl = [](auto&&self, std::size_t const n, std::size_t const k) -> std::size_t {
        switch (k) {
        case 0u:
          return 1u;
        case 1u:
          return n;
        case 2u:
          return (n * (n - 1u)) >> 1u;
        };
        return (self(self, n - 1u, k - 1u) * n) / k;
      };
      return impl(impl, n, k);
    };
    return (k > n) ? 0u : bin(n, std::min(k,n-k));
  }

  /**
   *  @brief Calculate binomial coeffient \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
   *         \f$\binom{n}{k}=\begin{cases}\frac{n^{k}}{k!} & \mathrm{if}\:k\leq\frac{n}{2}\\\frac{n^{\left(n-k\right)}}{\left(n-k\right)!} & \mathrm{if}\:k>\frac{n}{2}\end{cases}\f$ or falling factorial
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @return 0 otherwise
   */
  [[nodiscard]] static constexpr std::size_t binomialFallingFactorial(std::size_t const n, std::size_t const k) noexcept {
      auto constexpr impl = [](std::size_t const n, std::size_t const k) -> std::size_t {
          std::size_t const n2k = fallingFactorial(n, k);
          std::size_t const kFact = factorial(k);

          return n2k / kFact;
      };

      return ( k > n ) ? 0u : impl(n, std::min(k, n-k));
  }

  /**
   *  @brief Calculate binomial coeffient \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as naive
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @return 0 otherwise
   */
  [[nodiscard]] static constexpr std::size_t binomialNaive(std::size_t const n, std::size_t k) noexcept {
      auto constexpr impl = [](std::size_t const n, std::size_t const k) -> std::size_t {
          std::size_t const nF = factorial( n );
          std::size_t const nkF = factorial( n - k );
          std::size_t const kF = factorial( k );

          return (nF / nkF) / kF;
      };

      return ( k > n ) ? 0u : impl(n, std::min(k, n-k));
  }

  /**
   *  @brief Calculate binomial coeffient \f$\binom{n}{k}=\frac{n!}{k!\left(n-k\right)!}\f$ as
   *         \f$\binom{n}{k}=\binom{n-1}{k-1}+\binom{n-1}{k}\f$ with in place vector
   *
   *  @param n "k objects can be chosen from among n objects", n >= k >= 0
   *  @param k "k objects can be chosen from among n objects", n >= k >= 0
   *  @return \f$\binom{n}{k}\f$, if n >= k >= 0
   *  @return 0 otherwise
   */
  [[nodiscard]] static constexpr std::size_t binomialInPlace(std::size_t const n, std::size_t k) noexcept {
      if( k > n ) {
          return 0;
      }

      k = std::min(k,n-k);
      if( (n == 0) || (k==0)){
          return 1;
      }
      else if ( k==1){
          return n;
      }
      else{
          std::vector<std::size_t> column(n<<2u,0);
          std::size_t* prev = column.data();
          std::size_t* curr = prev + n;

          for( size_t i = 0; i <= n; ++i ){
              prev[i]=1;
          }

          for( std::size_t i_k = 1; i_k <= k; ++i_k) {
              curr[i_k]=1;
              for( std::size_t j_n = i_k+1; j_n < n; ++j_n ) {
                  curr[j_n] = curr[j_n - 1] + prev[j_n - 1];
              }

              std::swap(curr, prev);
          }

          return curr[n-1] + prev[n-1];
      }
  }

  /**
   *  @brief Return 'x' to integer power 'i' \f$x^{i},x\in\mathbb{R},i\in\mathbb{N}_{0}\f$
   *
   *  @param x real value
   *  @param i ith integer power
   *  @return x to ith power
   */
  [[nodiscard]] static constexpr real toIntegerPower(real x, std::size_t i) noexcept {
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

  /**
   *  @brief Return factorial n!
   *
   *  @param n integer value
   *  @return Factorial of n
   */
   [[nodiscard]] static constexpr std::size_t factorial(std::size_t const n) noexcept {
       switch(n){
       case 0:
       case 1:
           return 1;

       case 2:
           return 2;

       default:
           std::size_t product = 2;
           for( std::size_t i = 3; i <= n; ++i){
               product *= i;
           }
           return product;
        }
    }

  /**
   *  @brief Return falling factorial \f$n^{\underline{k}}=n\cdot\left(n-1\right)\cdot\ldots\cdot\left(n-\left(k-1\right)\right)\f$
   *
   *  @param n integer value, n >= k >= 0
   *  @param k integer value, n >= k >= 0
   *  @return falling factorial
   */
    [[nodiscard]] static constexpr std::size_t fallingFactorial(std::size_t const n, std::size_t const k) noexcept {
        if( ( n == 0 ) || ( k == 0 ) ) {
            return 1;
        } else {
            std::size_t prod = (n - k) + 1;
            for( size_t i = prod +1; i <= n; ++i ){
                prod *= i;
            }
            return prod;
        }
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
  [[nodiscard]] static constexpr real B(int const i, int const n, real const u) noexcept {
    if ((i > n) || (n < 0) || (i < 0)) {
      return real(0);
    }
    real const b = real(binomialRecursive(n, i)) * toIntegerPower(u, i) * toIntegerPower(real(1) - u, n - i);
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
  [[nodiscard]] static constexpr real dB(int const i, int const n, real const u) noexcept {
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
  [[nodiscard]] static constexpr real d2B(int const i, int const n, real const u) noexcept {
    real const b = real(n) * real(n - 1) * (B(i, n - 2, u) - real(2) * B(i - 1, n - 2, u) + B(i - 2, n - 2, u));
    return b;
  }
};

} // namespace curve::bezier::utilities

#endif
