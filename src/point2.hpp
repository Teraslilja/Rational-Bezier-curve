//
// (C) Matti Lehtonen 2023
//

/**
 *  @file point.hpp This file contains two implementations of point classes; Point2 and Point3 for 2D and 3D curves.
 */

#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <iostream>
#include <optional>
#include <type_traits>

namespace curve {

/**
 *  @brief A 2D point to be used with rational bezier curve
 */
template <typename type = float>
requires std::is_floating_point_v<type>
class Point2 {
public:
  using real = type;
  using point = Point2<real>;

  /**
   *  @brief Constructor with coordinates
   *
   *  @param x the X coordinate of point
   *  @param y the Y coordinate of point
   */
  inline constexpr Point2(real const x, real const y) noexcept : x_(x), y_(y) {}

  /**
   *  @brief The default constructor. Coordinated are zeroed
   */
  inline constexpr Point2() noexcept : Point2(real(0), real(0)) {}

  /**
   *  @brief The destructor
   */
  inline constexpr ~Point2() noexcept = default;

  /**
   *  @brief Return value of X coordinate
   *
   *  @return X coordinate for reading
   */
  [[nodiscard]] inline constexpr real x() const noexcept { return this->x_; }

  /**
   *  @brief Return reference to X coordinate for modifying purposes
   *
   *  @return reference to X coordinate for writing
   */
  [[nodiscard]] inline constexpr real &x() noexcept { return this->x_; }

  /**
   *  @brief Return value of Y coordinate
   *
   *  @return Y coordinate for reading
   */
  [[nodiscard]] inline constexpr real y() const noexcept { return this->y_; }

  /**
   *  @brief Return reference to Y coordinate for modifying purposes
   *
   *  @return reference to Y coordinate for writing
   */
  [[nodiscard]] inline constexpr real &y() noexcept { return this->y_; }

  /**
   *  @brief Multiply coordinate values of point with 'scale'
   *
   *  @param scale multiplication factor
   *  @return new point that is multiplied
   */
  inline constexpr point operator*(real const scale) const noexcept {
    return point(scale * this->x_, scale * this->y_);
  }

  /**
   *  @brief Dot product of two points
   *
   *  @param p another point
   *  @return Dot product of two points
   */
  inline constexpr real operator*(point const p) const noexcept { return this->x_ * p.x_ + this->y_ * p.y_; }

  /**
   *  @brief Divide coodinate vales of point by 'div'
   *
   *  @param div divider
   *  @return the divided point as std::optional
   *  @return std::nullopt, if divided by zero or near zero
   */
  [[nodiscard]] inline constexpr std::optional<point> operator/(real const div) const noexcept {
    bool const state = std::abs(div) < std::numeric_limits<real>::epsilon();
    return state ? std::nullopt : std::make_optional(point(this->x_ / div, this->y_ / div));
  }

  /**
   *  @brief Calculate substraction between this and another point 'p'
   *
   *  @param p point to be subracted from 'this'
   *  @return a new point containing the substraction
   */
  [[nodiscard]] inline constexpr point operator-(point p) const noexcept {
    return point(this->x() - p.x(), this->y() - p.y());
  }

  /**
   *  @brief Calculate addition between this and point 'p'
   *
   *  @param p point to be added to 'this'
   *  @return a new point containing the addition
   */
  [[nodiscard]] inline constexpr point operator+(point const p) const noexcept {
    return Point2(this->x_ + p.x_, this->y_ + p.y_);
  }

  /**
   *  @brief Add point 'p' to this point
   *
   *  @param p point to be added to 'this'
   *  @return 'this' point
   */
  inline constexpr point &operator+=(point const p) noexcept {
    this->x_ += p.x_;
    this->y_ += p.y_;
    return *this;
  }

  /**
   *  @brief Calculate trace of point or sum of coordinate values. This is different than 1-norm:
   * \f$\stackrel[i=1]{n}{\sum}\left|x_{i}\right|\f$
   *
   *  @return sum of coordinates
   */
  [[nodiscard]] inline constexpr real trace() const noexcept { return this->x_ + this->y_; }

  /**
   *  @brief Calculate squared length of point as vector
   *
   *  @return squared length
   */
  [[nodiscard]] inline constexpr real lengthSquared() const noexcept {
    return this->x_ * this->x_ + this->y_ * this->y_;
  }

  /**
   *  @brief Calcualate length of point as vector
   *
   *  @return length
   */
  [[nodiscard]] inline constexpr real length() const noexcept { return std::sqrt(this->lengthSquared()); }

  /**
   *  @brief Calculate distance between this point and point 'p'
   *
   *  @param p the point calcualte distance to
   *  @return distance between points
   */
  [[nodiscard]] inline constexpr real distance(point p) const noexcept {
    point const difference = *this - p;
    return difference.length();
  }

  /**
   *  @brief Normalize lenth of point (as vector) to unit length
   *
   *  @return a point, which length is one as std::optional
   *  @return std::nullopt, if length of point is zero or near zero
   */
  [[nodiscard]] inline constexpr std::optional<point> normalize() const noexcept {
    real const length = std::sqrt(this->x_ * this->x_ + this->y_ * this->y_);
    return *this / length;
  }

  /**
   *  @brief Return order or number of coordinate values point have
   *
   *  @return 2
   */
  [[nodiscard]] static inline consteval std::size_t getOrder() noexcept { return 2u; }

  /**
   *  @brief Output point 'data' to stream 'out'
   *
   *  @param out the output stream
   *  @param data the point to be streamed
   *  @return the 'out' stream
   */
  inline friend std::ostream &operator<<(std::ostream &out, point const &data) {
    out << "{";
    out << data.x_ << "," << data.y_;
    out << "}";
    return out;
  }

  /**
   *  @brief Output optional point 'data' to stream 'out'
   *
   *  @param out the output stream
   *  @param data the point to be streamed, or if std::nullopt, string "'no value'"
   *  @return the 'out' stream
   */
  inline friend std::ostream &operator<<(std::ostream &out, std::optional<point> const &data) {
    if (data.has_value()) {
      out << data.value();
    } else {
      out << "'no value'";
    }
    return out;
  }

private:
  //! The coordinate values @{
  real x_, y_;
  //! @}
};

/**
 *  @brief A 3D point to be used with rational bezier curve
 */
template <typename type = float>
requires std::is_floating_point_v<type>
class Point3 {
public:
  using real = type;
  using point = Point3<real>;

  /**
   *  @brief Constructor with coordinates
   *
   *  @param x the X coordinate of point
   *  @param y the Y coordinate of point
   *  @param z the Z coordinate of point
   */
  inline constexpr Point3(real const x, real const y, real const z) noexcept : x_(x), y_(y), z_(z) {}

  /**
   *  @brief The default constructor. Coordinated are zeroed
   */
  inline constexpr Point3() noexcept : Point3(real(0), real(0), real(0)) {}

  /**
   *  @brief Construct 3D point from 2D point and zero Z coordinate
   *
   *  @param p 2D point
   */
  inline constexpr Point3(Point2<real> const p) noexcept : Point3(p.x(), p.y(), real(0)) {}

  /**
   *  @brief The destructor
   */
  inline constexpr ~Point3() noexcept = default;

  /**
   *  @brief Return value of X coordinate
   *
   *  @return X coordinate for reading
   */
  [[nodiscard]] inline constexpr real x() const noexcept { return this->x_; }

  /**
   *  @brief Return reference to X coordinate for modifying purposes
   *
   *  @return reference to X coordinate for writing
   */
  [[nodiscard]] inline constexpr real &x() noexcept { return this->x_; }

  /**
   *  @brief Return value of Y coordinate
   *
   *  @return Y coordinate for reading
   */
  [[nodiscard]] inline constexpr real y() const noexcept { return this->y_; }

  /**
   *  @brief Return reference to Y coordinate for modifying purposes
   *
   *  @return reference to Y coordinate for writing
   */
  [[nodiscard]] inline constexpr real &y() noexcept { return this->y_; }

  /**
   *  @brief Return value of Z coordinate
   *
   *  @return Z coordinate for reading
   */
  [[nodiscard]] inline constexpr real z() const noexcept { return this->z_; }

  /**
   *  @brief Return reference to Z coordinate for modifying purposes
   *
   *  @return reference to Z coordinate for writing
   */
  [[nodiscard]] inline constexpr real &z() noexcept { return this->z_; }

  /**
   *  @brief Multiply coordinate values of point with 'scale'
   *
   *  @param scale multiplication factor
   *  @return new point that is multiplied
   */
  [[nodiscard]] inline constexpr point operator*(real const scale) const noexcept {
    return point(scale * this->x_, scale * this->y_, scale * this->z_);
  }

  /**
   *  @brief Dot product of two points
   *
   *  @param p multiplication factor
   *  @return dot product of two points
   */
  inline constexpr real operator*(Point2<real> const p) const noexcept {
    return this->x_ * p.x_ + this->y_ * p.y_ + this->z_ * p.z_;
  }

  /**
   *  @brief Divide coodinate vales of point by 'div'
   *
   *  @param div divider
   *  @return the divided point as std::optional
   *  @return std::nullopt, if divided by zero or near zero
   */
  [[nodiscard]] inline constexpr std::optional<point> operator/(real const div) const noexcept {
    bool const state = std::abs(div) < std::numeric_limits<real>::epsilon();
    return state ? std::nullopt : std::make_optional(point(this->x_ / div, this->y_ / div, this->z_ / div));
  }

  /**
   *  @brief Calculate substraction between this and another point 'p'
   *
   *  @param p point to be subracted from 'this'
   *  @return a new point containing the substraction
   */
  [[nodiscard]] inline constexpr point operator-(point p) const noexcept {
    return point(this->x_ - p.x_, this->y_ - p.y_, this->z_ - p.z_);
  }

  /**
   *  @brief Calculate addition between this and point 'p'
   *
   *  @param p point to be added to 'this'
   *  @return a new point containing the addition
   */
  [[nodiscard]] inline constexpr point operator+(point const p) const noexcept {
    return Point3(this->x_ + p.x_, this->y_ + p.y_, this->z_ + p.z_);
  }

  /**
   *  @brief Add point 'p' to this point
   *
   *  @param p point to be added to 'this'
   *  @return 'this' point
   */
  inline constexpr point &operator+=(point const p) noexcept {
    this->x_ += p.x_;
    this->y_ += p.y_;
    this->z_ += p.z_;
    return *this;
  }

  /**
   *  @brief Calculate trace of point or sum of coordinate values. This is
   * different than 1-norm: \f$\stackrel[i=1]{n}{\sum}\left|x_{i}\right|\f$
   *
   *  @return sum of coordinates
   */
  [[nodiscard]] inline constexpr real trace() const noexcept { return this->x_ + this->y_ + this->z_; }

  /**
   *  @brief Calculate squared length of point as vector
   *
   *  @return squared length
   */
  [[nodiscard]] inline constexpr real lengthSquared() const noexcept {
    return this->x_ * this->x_ + this->y_ * this->y_ + this->z_ * this->z_;
  }

  /**
   *  @brief Calcualate length of point as vector
   *
   *  @return length
   */
  [[nodiscard]] inline constexpr real length() const noexcept { return std::sqrt(this->lengthSquared()); }

  /**
   *  @brief Calculate distance between this point and point 'p'
   *
   *  @param p the point calcualte distance to
   *  @return distance between points
   */
  [[nodiscard]] inline constexpr real distance(point p) const noexcept {
    point const difference = *this - p;
    return difference.length();
  }

  /**
   *  @brief Normalize lenth of point (as vector) to unit length
   *
   *  @return a point, which length is one as std::optional
   *  @return std::nullopt, if length of point is zero or near zero
   */
  [[nodiscard]] inline constexpr std::optional<point> normalize() const noexcept {
    real const length = std::sqrt(this->x_ * this->x_ + this->y_ * this->y_ + this->z_ * this->z_);
    return *this / length;
  }

  /**
   *  @brief Return order or number of coordinate values point have
   *
   *  @return 3
   */
  [[nodiscard]] static inline consteval std::size_t getOrder() noexcept { return 3u; }

  /**
   *  @brief Output point 'data' to stream 'out'
   *
   *  @param out the output stream
   *  @param data the point to be streamed
   *  @return the 'out' stream
   */
  inline friend std::ostream &operator<<(std::ostream &out, point const &data) {
    out << "{";
    out << data.x_ << "," << data.y_ << "," << data.z_;
    out << "}";
    return out;
  }

  /**
   *  @brief Output optional point 'data' to stream 'out'
   *
   *  @param out the output stream
   *  @param data the point to be streamed, or if std::nullopt, string "'no value'"
   *  @return the 'out' stream
   */
  inline friend std::ostream &operator<<(std::ostream &out, std::optional<point> const &data) {
    if (data.has_value()) {
      out << data.value();
    } else {
      out << "'no value'";
    }
    return out;
  }

private:
  //! The coordinate values @{
  real x_, y_, z_;
  //! @}
};

} // namespace curve

#endif
