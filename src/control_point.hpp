//
// (C) Matti Lehtonen 2023
//

/**
 *  @file control_point.hpp This file contains module interface for control point
 */

#ifndef CONTROL_POINT_H
#define CONTROL_POINT_H

#include "point2.hpp"
#include "point3.hpp"

#include <vector>

namespace curve::bezier::rational {

/**
 *  @brief Control point class for rational bezier curves. Uses either @ref curve::points::Point2 or @ref
 * curve::points::Point3
 */
template <class P, typename type = typename P::real>
requires(std::is_same_v<P, curve::points::Point2<typename P::real>> ||
         std::is_same_v<P, curve::points::Point3<typename P::real>>) &&
    std::is_floating_point_v<type> class ControlPoint {
public:
  using Point = P;
  using Controlpoint = ControlPoint<Point, typename Point::real>;
  using real = type;

  /**
   *  @brief Control point constructor with weight and point
   *
   *  @param w the weight of point
   *  @param p the point
   */
  inline constexpr ControlPoint(real const w, Point const p) noexcept : w_(w), p_(p) {}

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
  [[nodiscard]] inline constexpr Point p() const noexcept { return this->p_; }

  /**
   *  @brief Return reference to (coordinate) values of control point for modifying purposes
   *
   *  @return the reference to point
   */
  [[nodiscard]] inline constexpr Point &p() noexcept { return this->p_; }

  /**
   *  @brief Stream control point 'data' to output stream 'out'
   *
   *  @param out output stream
   *  @param data control point
   *  @return 'out' stream
   */
  friend std::ostream &operator<<(std::ostream &out, Controlpoint const &data) {
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
  friend std::ostream &operator<<(std::ostream &out, std::optional<Controlpoint> const &data) {
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
  friend std::ostream &operator<<(std::ostream &out, std::vector<Controlpoint> const &data) {
    out << "[";
    for (std::size_t i = 0u; auto const &p : data) {
      out << ((i > 0u) ? ", " : "") << p;

      ++i;
    }
    out << "]";
    return out;
  }

protected:
  real w_;  //< Weight of control point of rational bezier curve
  Point p_; //< Point of control point of rational bezier curve
};

} // namespace curve::bezier::rational

#endif
