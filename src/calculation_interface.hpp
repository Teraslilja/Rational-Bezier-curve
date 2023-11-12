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
namespace {

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
   *  @brief Store a pair of location of curve and a point
   */
  struct PointAtCurve {
    real const u;     //< Location at curve
    ConstPoint point; //< Point at curve

    /**
     *  @brief Constructor
     */
    inline constexpr PointAtCurve(real const u, ConstPoint point) noexcept : u{u}, point{point} {}
  };

  /**
   *  @brief Store a pair of squared distance and location of curve
   */
  struct DistanceFromLocation {
    real distanceSquared; //< Squared distance from the curve at location u
    real u;               //< Location at curve

    /// @brief Rejected distance
    static constexpr real const REJECTED_DISTANCE = std::numeric_limits<real>::infinity();

    /// @brief Value marked to indicate invalid value of u
    static constexpr real const INVALID_U = std::numeric_limits<real>::infinity();

    /// @brief Constructor
    inline constexpr DistanceFromLocation(real const distanceSquared, real const u) noexcept
        : distanceSquared{distanceSquared}, u{u} {}

    /// @brief A default constructor
    inline constexpr DistanceFromLocation() noexcept : DistanceFromLocation(REJECTED_DISTANCE, INVALID_U) {}

    /// @brief A copy constructor
    inline constexpr DistanceFromLocation(DistanceFromLocation const &v) noexcept
        : DistanceFromLocation(v.distanceSquared, v.u) {}

    /// @brief A copy operator
    inline constexpr DistanceFromLocation &operator=(DistanceFromLocation const &v) noexcept {
      this->distanceSquared = v.distanceSquared;
      this->u = v.u;

      return *this;
    };

    /**
     *  @brief Is the location of curve a valid one
     *
     *  @return true, if u is not infinity
     *  @return false, if u is infinity
     */
    inline constexpr bool hasValidLocation() const noexcept { return std::isfinite(this->u); }
  };

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
  constexpr void approximateAsLinestring(std::size_t const initial_vertices, real const max_segment_error,
                                         std::function<void(real const, PointAtCurve const)> pick_segment) const;
  constexpr void initialGuessesFromCurve(std::vector<DistanceFromLocation> &nearest, ConstPoint p,
                                         real const curveApproximation) const noexcept;
  [[nodiscard]] constexpr Point findNearestPointFor(ConstPoint p) const;
  /// @}
};

} // namespace
} // namespace rational
} // namespace bezier
} // namespace curve

#endif
