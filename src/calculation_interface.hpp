//
// (C) Matti Lehtonen 2023
//

/**
 *  @file calculation_interface.hpp This file contains module interface for calculation interface
 */

#ifndef CALCULATION_INTERFACE_H
#define CALCULATION_INTERFACE_H

#include "control_point.hpp"

#include <span>

namespace curve::bezier::rational {

/**
 *  @brief Define the calculation interface
 */
template <class CP>
requires std::is_same_v<CP, curve::bezier::rational::ControlPoint<typename CP::Point, typename CP::real>>
struct CalculationInterface {
  using ControlPoint = CP;
  using Point = typename ControlPoint::Point;
  using real = typename Point::real;

  using ControlPointSpan = typename std::span<ControlPoint const>;

  /**
   *  @brief A default constructor
   */
  constexpr CalculationInterface() = default;

  /**
   *  @brief A default destructor
   */
  constexpr ~CalculationInterface() = default;

  /// @brief The methods of calculation interface @{
  [[nodiscard]] constexpr std::size_t numberOfControlPoints() const noexcept;
  [[nodiscard]] constexpr ControlPointSpan const &getSpan() const noexcept;

  [[nodiscard]] constexpr Point C(real const u) const noexcept;
  [[nodiscard]] constexpr Point dC(real const u) const noexcept;
  [[nodiscard]] constexpr Point d2C(ControlPointSpan const controlPoints, real const u) noexcept;

  [[nodiscard]] constexpr real curveLength() const noexcept;

  [[nodiscard]] constexpr std::vector<Point> asLinestring() const;

  [[nodiscard]] constexpr Point findNearestPointFor(Point const p) const;
  /// @}
};

} // namespace curve::bezier::rational

#endif
