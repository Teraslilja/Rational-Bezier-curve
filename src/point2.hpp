//
// (C) Matti Lehtonen 2023
//

/**
 *  @file point2.hpp This file contains module interface for 2D point
 */

#ifndef POINT2_H
#define POINT2_H

#include "point_base.hpp"

#include <cmath>
#include <ostream>

namespace curve::points {

/**
 *  @brief A 2D point to be used with curves
 */
template<typename type = float>
  requires std::is_floating_point_v<type>
class Point2 : private PointD<2u, type>
{
private:
  using Base = PointD<2u, type>;

protected:
  /// @brief Define the dimensions
  enum class Dimension : std::size_t
  {
    X,
    Y
  };

public:
  using real = type;
  using Point = Point2<real>;

public:
  /**
   *  @brief The default constructor. Coordinated are zeroed
   */
  constexpr Point2() noexcept
    : Point2(real(0), real(0))
  {
  }

  /**
   *  @brief The destructor
   */
  constexpr ~Point2() noexcept = default;

  /**
   *  @brief Constructor with coordinates
   *
   *  @param x the x coordinates of point
   *  @param y the y coordinates of point
   */
  constexpr Point2(real const x, real const y) noexcept
  {
    this->x() = x;
    this->y() = y;
  }

public:
  /**
   *  @brief Return value of X coordinate
   *
   *  @return X coordinate for reading
   */
  [[nodiscard]] constexpr real x() const noexcept
  {
    constexpr const std::size_t index = std::size_t(Dimension::X);
    static_assert(index < this->getOrder());
    return this->coords_[index];
  }

  /**
   *  @brief Return reference to X coordinate for modifying purposes
   *
   *  @return reference to X coordinate for writing
   */
  [[nodiscard]] constexpr real& x() noexcept
  {
    constexpr const std::size_t index = std::size_t(Dimension::X);
    static_assert(index < this->getOrder());
    return this->coords_[index];
  }

  /**
   *  @brief Return value of Y coordinate
   *
   *  @return Y coordinate for reading
   */
  [[nodiscard]] constexpr real y() const noexcept
  {
    constexpr const std::size_t index = std::size_t(Dimension::Y);
    static_assert(index < this->getOrder());
    return this->coords_[index];
  }

  /**
   *  @brief Return reference to Y coordinate for modifying purposes
   *
   *  @return reference to Y coordinate for writing
   */
  [[nodiscard]] constexpr real& y() noexcept
  {
    constexpr const std::size_t index = std::size_t(Dimension::Y);
    static_assert(index < this->getOrder());
    return this->coords_[index];
  }

  /**
   *  @brief Multiply coordinate values of point with 'scale'
   *
   *  @param scale multiplication factor
   *  @return new point that is multiplied
   */
  constexpr Point operator*(real const scale) const noexcept
  {
    return Point(scale * this->x(), scale * this->y());
  }

  /**
   *  @brief Dot product of two points
   *
   *  @param p another point
   *  @return Dot product of two points
   */
  constexpr real operator*(Point const p) const noexcept
  {
    return this->x() * p.x() + this->y() * p.y();
  }

  /**
   *  @brief Divide coodinate vales of point by 'div'
   *
   *  @param div divider
   *  @return the divided point as std::optional
   *  @return std::nullopt, if divided by zero or near zero
   */
  [[nodiscard]] constexpr std::optional<Point> operator/(
    real const div) const noexcept
  {
    bool const state = std::abs(div) < std::numeric_limits<real>::epsilon();
    return state ? std::nullopt
                 : std::make_optional(Point(this->x() / div, this->y() / div));
  }

  /**
   *  @brief Calculate substraction between this and another point 'p'
   *
   *  @param p point to be subracted from 'this'
   *  @return a new point containing the substraction
   */
  [[nodiscard]] constexpr Point operator-(Point const p) const noexcept
  {
    return Point(this->x() - p.x(), this->y() - p.y());
  }

  /**
   *  @brief Calculate addition between this and point 'p'
   *
   *  @param p point to be added to 'this'
   *  @return a new point containing the addition
   */
  [[nodiscard]] constexpr Point operator+(Point const p) const noexcept
  {
    return Point(this->x() + p.x(), this->y() + p.y());
  }

  /**
   *  @brief Add point 'p' to this point
   *
   *  @param p point to be added to 'this'
   *  @return 'this' point
   */
  constexpr Point& operator+=(Point const p) noexcept
  {
    this->x() += p.x();
    this->y() += p.y();
    return *this;
  }

  /**
   *  @brief Calculate trace of point or sum of coordinate values. This is
   * different than 1-norm: \f$\stackrel[i=1]{n}{\sum}\left|x(){i}\right|\f$
   *
   *  @return sum of coordinates
   */
  [[nodiscard]] constexpr real trace() const noexcept
  {
    return this->x() + this->y();
  }

  /**
   *  @brief Calculate squared length of point as vector
   *
   *  @return squared length
   */
  [[nodiscard]] constexpr real lengthSquared() const noexcept
  {
    return this->x() * this->x() + this->y() * this->y();
  }

  /**
   *  @brief Calcualate length of point as vector
   *
   *  @return length
   */
  [[nodiscard]] constexpr real length() const noexcept
  {
    return std::sqrt(this->lengthSquared());
  }

  /**
   *  @brief Calculate distance between this point and point 'p'
   *
   *  @param p the point calcualte distance to
   *  @return distance between points
   */
  [[nodiscard]] constexpr real distance(Point const p) const noexcept
  {
    Point const difference = *this - p;
    return difference.length();
  }

  /**
   *  @brief Normalize lenth of point (as a vector) to unit length
   *
   *  @return a point, which length is one as std::optional
   *  @return std::nullopt, if length of point is zero or near zero
   */
  [[nodiscard]] constexpr std::optional<Point> normalize() const noexcept
  {
    return *this / this->length();
  }

  /**
   *  @brief Return order or number of coordinate values point have
   *
   *  @return DIM
   */
  [[nodiscard]] static inline consteval std::size_t getOrder() noexcept
  {
    return Base::dim;
  }

  /**
   *  @brief Stream point 'data' to output stream 'out'
   *
   *  @param out the output stream
   *  @param data the point to be streamed
   *  @return the 'out' stream
   */
  inline friend std::ostream& operator<<(std::ostream& out, Point const& data)
  {
    out << "{";
    out << data.x() << "," << data.y();
    out << "}";
    return out;
  }

  /**
   *  @brief Stream optional point 'data' to output stream 'out'
   *
   *  @param out the output stream
   *  @param data the point to be streamed, or if std::nullopt, string "'no
   * value'"
   *  @return the 'out' stream
   */
  inline friend std::ostream& operator<<(std::ostream& out,
                                         std::optional<Point> const& data)
  {
    if (data.has_value()) {
      out << data.value();
    } else {
      out << "'no value'";
    }
    return out;
  }
};

} // namespace curve::points

#endif
