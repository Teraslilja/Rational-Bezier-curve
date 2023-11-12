//
// (C) Matti Lehtonen 2023
//

/**
 *  @file control_point.hpp This file contains implementation of control point class.
 */

#ifndef CONTROL_POINT_H
#define CONTROL_POINT_H

#include "point.hpp"

#include <vector>

namespace curve {
namespace bezier {
namespace rational {

/**
 *  @brief Control point class for rational bezier curves. Uses either @ref curve::Point2 or @ref curve::Point3
 */
template <class P, typename type = typename P::real>
requires(std::is_same_v<P, Point2<typename P::real>> || std::is_same_v<P, Point3<typename P::real>>) &&
    std::is_floating_point_v<type> class ControlPoint {
public:
  using point = P;
  using real = type;

  /**
   *  @brief Control point constructor with weight and point
   *
   *  @param w the weight of point
   *  @param p the point
   */
  inline constexpr ControlPoint(real const w, point const p) noexcept : w_(w), p_(p) {}

  /**
   *  @brief Return value of weight for reading purposes
   *
   *  @return the weight of control point
   */
  [[nodiscard]] inline constexpr real w() const noexcept { return this->w_; }

  /**
   *  @brief Return reference of value of weight for modifying purposes
   *
   *  @return a reference to weight of control point
   */
  [[nodiscard]] inline constexpr real &w() noexcept { return this->w_; }

  /**
   *  @brief Return (coordinate) values of control point for reading purposes
   *
   *  @return the point
   */
  [[nodiscard]] inline constexpr point p() const noexcept { return this->p_; }

  /**
   *  @brief Return reference to (coordinate) values of control point for modifying purposes
   *
   *  @return the reference to point
   */
  [[nodiscard]] inline constexpr point &p() noexcept { return this->p_; }

  /**
   *  @brief Stream control point 'data' to output stream 'out'
   *
   *  @param out output stream
   *  @param data control point
   *  @return 'out' stream
   */
  friend inline std::ostream &operator<<(std::ostream &out, ControlPoint<point> const &data) {
    out << "{";
    out << data.w_ << "," << data.p_;
    out << "}";
    return out;
  }

  /**
   *  @brief Stream optional control point 'data' to output stream 'out'
   *
   *  @param out output stream
   *  @param data optional control point or if std::nullopt output string "'no value'"
   *  @return 'out' stream
   */
  friend inline std::ostream &operator<<(std::ostream &out, std::optional<ControlPoint<point>> const &data) {
    if (data.has_value()) {
      out << data.value();
    } else {
      out << "'no value'";
    }
    return out;
  }

  /**
   *  @brief Stream a vector of control points 'data' to output stream 'out'
   *
   *  @param out output stream
   *  @param data vector of control points
   *  @return 'out' stream
   */
  friend inline std::ostream &operator<<(std::ostream &out, std::vector<ControlPoint> const &data) {
    out << "[";
    if (!data.empty()) {
      out << data.front();
      for (auto p = std::next(data.cbegin()); p != data.cend(); ++p) {
        out << ", " << *p;
      }
    }
    out << "]";
    return out;
  }

protected:
  real w_;  //< Weight of control point of rational bezier curve
  point p_; //< Point of control point of rational bezier curve
};

} // namespace rational
} // namespace bezier
} // namespace curve
#endif
