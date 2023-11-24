//
// (C) Matti Lehtonen 2023
//

/**
 *  @file validation_interface.hpp This file contains module interface for validation interface
 */

#ifndef VALIDATION_INTERFACE_H
#define VALIDATION_INTERFACE_H

#include "control_point.hpp"

#include <cstdint>
#include <variant>

namespace curve::bezier::rational {

/**
 *  @brief Isses that rational bezier curve may have
 */
enum class ValidityIssue : std::uint16_t {
  ISSUE_NOT_ENOUGHT_CONTROL_POINTS = 1u, //< Not enought control points to calculate
  ISSUE_U_IS_INVALID,                    //< Parameter U is not in valid range [0;1]
  ISSUE_BAD_COMBINATION_OF_WEIGHTS,      //< Sum of weights is (too near) zero
  ISSUE_OUT_OF_HEAP_MEMORY,              //< Container (re)allocation has failed
  ISSUE_BAD_CONTROLPOINT_WEIGHT,         //< Weight of control point is infinite or NaN
  ISSUE_BAD_POINT,                       //< Coordiante value of point is infinite or NaN
};

/**
 *  @brief Output enumerated issue number 'data' to stream 'out'
 *
 *  @param out the output stream
 *  @param data the issue to be streamed
 *  @return the 'out' stream
 */
std::ostream &operator<<(std::ostream &out, ValidityIssue const data);

/**
 *  @brief Output optional enumerated issue number 'data' to stream 'out'
 *
 *  @param out the output stream
 *  @param data the issue to be streamed or "no value" string
 *  @return the 'out' stream
 *
 */
std::ostream &operator<<(std::ostream &out, std::optional<ValidityIssue> const &data);

/**
 *  @brief The delegation interface between classes ValidateRational and Rational
 */
template <class CP>
requires std::is_same_v<CP, ControlPoint<typename CP::Point, typename CP::real>>
struct ValidationInterface {
public:
  using ControlPoint = CP;
  using Point = typename ControlPoint::Point;
  using real = typename Point::real;

  using point_or_issue = std::variant<ValidityIssue, Point>;
  using real_or_issue = std::variant<ValidityIssue, real>;
  using vector_of_points_or_issue = std::variant<ValidityIssue, std::vector<Point>>;

public:
  /**
   *  @brief A default constructor
   */
  inline constexpr ValidationInterface() = default;

  /**
   *  @brief A default destructor
   */
  inline constexpr ~ValidationInterface() = default;

  /// The methods of delegated interface @{
  [[nodiscard]] constexpr point_or_issue pointAt(real const u) const noexcept;
  [[nodiscard]] constexpr real_or_issue curveLength() const noexcept;
  [[nodiscard]] constexpr vector_of_points_or_issue asLineString() const noexcept;
  [[nodiscard]] constexpr point_or_issue velocityAt(real const u) const noexcept;
  [[nodiscard]] constexpr real_or_issue speedAt(real const u) const noexcept;
  [[nodiscard]] constexpr point_or_issue tangentAt(real const u) const noexcept;
  [[nodiscard]] constexpr point_or_issue closestCurvePointFor(Point const p) const noexcept;
  /// @}
};

} // namespace curve::bezier::rational

#endif
