//
// (C) Matti Lehtonen 2023
//

/**
 *  @file newton_raphson.hpp This file contains module interface for
 * Newton-Raphson algorithm
 */

#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include <functional>
#include <optional>
#include <type_traits>

#include <cmath>

namespace curve::bezier::utilities {

/**
 *  @brief A structure to find root of function f(x) with help of f'(x)
 */
template<typename type = float>
  requires std::is_floating_point_v<type>
struct NewtonRaphson
{
public:
  using real = type;

public:
  /**
   *  @brief Find root of function f, f(x) = 0 by Newton-Raphson method.
   *  Iterate: \f$x_{i+1}=x_{i}-\frac{f(x_{i})}{f'(x_{i})}\f$
   *
   *  @param lowerBound if updated coordiante goes below this, reject iteration
   *  @param x initial starting value
   *  @param upperBound if updated coordiante goes aboce this, reject iteration
   *  @param roundLimit the maximum number of iteration rounds
   *  @param f function or lambda function to calculate f(x)
   *  @param df function or lambda function to calculate f'(x)
   *  @return std::nullopt, if cannot calculate f'(x)
   *  @return std::nullopt, if cannot calculate f(x)
   *  @return std::nullopt, new x_(i+1) do not belong range [ lowerBound ;
   * upperBound ]
   *  @return x, if |f'(x)| < DF_ZERO_EPS
   *  @return x otherwise as indication for f(x)=0
   */
  [[nodiscard]] static constexpr std::optional<real> findRoot(
    real const lowerBound,
    real x,
    real const upperBound,
    std::size_t const roundLimit,
    std::function<std::optional<real>(real const)> const f,
    std::function<std::optional<real>(real const)> const df)
  {
    real constexpr DF_ZERO_EPS = real(1e-5);
    real constexpr DX_ZERO_EPS = real(1e-6);

    for (std::size_t i = 0u; i < roundLimit; ++i) {
      std::optional<real> const df_x = df(x);
      if (!df_x.has_value()) {
        return std::nullopt;
      }

      std::optional<real> const f_x = f(x);
      if (!f_x.has_value()) {
        return std::nullopt;
      }

      if (std::abs(df_x.value()) < DF_ZERO_EPS) {
        // A local(?) minima/maxima or saddle point at f(u)
        return x;
      }

      real const delta_x = f_x.value() / df_x.value();
      if (std::abs(delta_x) < DX_ZERO_EPS) {
        return x;
      }

      // Update x for the next round
      x -= delta_x;
      if ((x < lowerBound) || (x > upperBound)) {
        // The method throws the search out of bounds
        return std::nullopt;
      }
    }

    return std::nullopt;
  }
};

} // namespace curve::bezier::utilities

#endif
