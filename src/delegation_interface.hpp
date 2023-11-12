//
// (C) Matti Lehtonen 2023
//

/**
 *  @file rational_bezier.hpp This file contains implementation of rational bezier class.
 */

#ifndef DELEGATION_INTERFACE_H
#define DELEGATION_INTERFACE_H

#include "control_point.hpp"

#include <cstdint>
#include <variant>

namespace curve {
namespace bezier {
namespace rational {

/**
 *  @brief Isses that rational bezier curve may have
 */
enum class ValidityIssue : std::uint16_t {
  ISSUE_NOT_ENOUGHT_CONTROL_POINTS = 1u, //< Not enought control points to calculate
  ISSUE_U_IS_INVALID,                    //< Parameter U is not in valid range [0;1]
  ISSUE_BAD_COMBINATION_OF_WEIGHTS,      //< Sum of weights is (too near) zero
  ISSUE_OUT_OF_HEAP_MEMORY,              //< Container (re)allocation has failed
};

// LCOV_EXCL_START
/**
 *  @brief Output enumerated issue number 'data' to stream 'out'
 *
 *  @param out the output stream
 *  @param data the issue to be streamed
 *  @return the 'out' stream
 */
inline std::ostream &operator<<(std::ostream &out, ValidityIssue const data) {
  switch (data) {
  case ValidityIssue::ISSUE_U_IS_INVALID:
    out << "ISSUE_U_IS_INVALID";
    break;
  case ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS:
    out << "ISSUE_NOT_ENOUGHT_CONTROL_POINTS";
    break;
  case ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS:
    out << "ISSUE_BAD_COMBINATION_OF_WEIGHTS";
    break;
  case ValidityIssue::ISSUE_OUT_OF_HEAP_MEMORY:
    out << "ISSUE_OUT_OF_HEAP_MEMORY";
    break;
  default:
    out << "Unknown issue (" << data << ")";
    break;
  }
  return out;
}

/**
 *  @brief Output optional enumerated issue number 'data' to stream 'out'
 *
 *  @param out the output stream
 *  @param data the issue to be streamed or "no value" string
 *  @return the 'out' stream
 *
 */
inline std::ostream &operator<<(std::ostream &out, std::optional<ValidityIssue> const &data) {
  if (data.has_value()) {
    out << data.value();
  } else {
    out << "'no value'";
  }
  return out;
}
// LCOV_EXCL_STOP

/**
 *  @brief The delegation interface between classes ValidateRational and Rational
 */
template <class CP>
requires std::is_same_v<CP, ControlPoint<typename CP::Point>>
struct DelegationInterface {
  using ControlPoint = CP;
  using ConstControlPoint = ControlPoint const;
  using Point = typename ControlPoint::Point;
  using ConstPoint = Point const;
  using real = typename Point::real;

  using point_or_issue = std::variant<ValidityIssue, Point>;
  using real_or_issue = std::variant<ValidityIssue, real>;
  using vector_of_points_or_issue = std::variant<ValidityIssue, std::vector<Point>>;

  /**
   *  @brief A default constructor
   */
  inline constexpr DelegationInterface() = default;

  /**
   *  @brief A default destructor
   */
  inline constexpr ~DelegationInterface() = default;

  /// The methods of delegated interface @{
  [[nodiscard]] constexpr point_or_issue pointAt(real const u) const noexcept;
  [[nodiscard]] constexpr real_or_issue curveLength() const noexcept;
  [[nodiscard]] constexpr vector_of_points_or_issue asLineString() const noexcept;
  [[nodiscard]] constexpr point_or_issue velocityAt(real const u) const noexcept;
  [[nodiscard]] constexpr real_or_issue speedAt(real const u) const noexcept;
  [[nodiscard]] constexpr point_or_issue tangentAt(real const u) const noexcept;
  [[nodiscard]] constexpr point_or_issue closestCurvePointFor(ConstPoint p) const noexcept;
  /// @}
};

} // namespace rational
} // namespace bezier
} // namespace curve

#endif
