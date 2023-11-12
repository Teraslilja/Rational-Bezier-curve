//
// (C) Matti Lehtonen 2023
//

/**
 *  @file rational_bezier.hpp This file contains implementation of rational bezier class.
 */

#ifndef CALCULATION_INTERFACE_H
#define CALCULATION_INTERFACE_H

#include "control_point.hpp"

#include <cstddef>
#include <functional>
#include <span>

namespace curve {
namespace bezier {
namespace rational {
namespace internal {

/**
 *  @brief Define the calculation interface
 */
template <class CP>
requires std::is_same_v<CP, ControlPoint<typename CP::point>>
struct CalculationInterface {
  using ControlPoint = CP; //< This class uses only read-only control points
  using point = typename ControlPoint::point;
  using real = typename point::real;
  using span = typename std::span<ControlPoint const>;

  /**
   *  @brief A default constructor
   */
  inline constexpr CalculationInterface() = default;

  /**
   *  @brief A default destructor
   */
  inline constexpr ~CalculationInterface() = default;

  /// The methods of delegated interface @{
  [[nodiscard]] constexpr point C(real const u) const noexcept;
  [[nodiscard]] inline constexpr std::size_t numberOfControlPoints() const noexcept;
  [[nodiscard]] inline constexpr span const &getSpan() const noexcept;
  [[nodiscard]] constexpr point dC(real const u) const noexcept;
  [[nodiscard]] constexpr point d2C(span const controlPoints, real const u) noexcept;
  [[nodiscard]] constexpr real length() const noexcept;
  [[nodiscard]] constexpr std::vector<point> asLinestring() const;
  constexpr void
  approximateAsLinestring(std::size_t const initial_vertices, real const max_segment_error,
                          std::function<void(real const, std::pair<real, point> const)> pick_segment) const;
  constexpr void initialGuessesFromCurve(std::vector<std::pair<real, real>> &nearest, point const p,
                                         real const curveApproximation) const noexcept;
  [[nodiscard]] constexpr point findNearestPointFor(point const p) const;
  /// @}
};

} // namespace internal
} // namespace rational
} // namespace bezier
} // namespace curve

#endif
