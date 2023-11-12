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
requires std::is_same_v<CP, ControlPoint<typename CP::Point>>
struct CalculationInterface {
  using ControlPoint = CP;
  using ConstControlPoint = ControlPoint const;
  using Point = typename ControlPoint::Point;
  using ConstPoint = Point const;
  using real = typename Point::real;

  using ControlPointSpan = typename std::span<ConstControlPoint>;

  /**
   *  @brief A default constructor
   */
  inline constexpr CalculationInterface() = default;

  /**
   *  @brief A default destructor
   */
  inline constexpr ~CalculationInterface() = default;

  /// The methods of delegated interface @{
  [[nodiscard]] constexpr Point C(real const u) const noexcept;
  [[nodiscard]] inline constexpr std::size_t numberOfControlPoints() const noexcept;
  [[nodiscard]] inline constexpr ControlPointSpan const &getSpan() const noexcept;
  [[nodiscard]] constexpr Point dC(real const u) const noexcept;
  [[nodiscard]] constexpr Point d2C(ControlPointSpan const controlPoints, real const u) noexcept;
  [[nodiscard]] constexpr real length() const noexcept;
  [[nodiscard]] constexpr std::vector<Point> asLinestring() const;
  constexpr void
  approximateAsLinestring(std::size_t const initial_vertices, real const max_segment_error,
                          std::function<void(real const, std::pair<real, ConstPoint> const)> pick_segment) const;
  constexpr void initialGuessesFromCurve(std::vector<std::pair<real, real>> &nearest, ConstPoint p,
                                         real const curveApproximation) const noexcept;
  [[nodiscard]] constexpr Point findNearestPointFor(ConstPoint p) const;
  /// @}
};

} // namespace internal
} // namespace rational
} // namespace bezier
} // namespace curve

#endif
