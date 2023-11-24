//
// (C) Matti Lehtonen 2023
//

/**
 *  @file point_base.hpp This file contains module interface for point base
 */

#ifndef POINT_BASE_H
#define POINT_BASE_H

#include <array>
#include <optional>
#include <type_traits>

namespace curve::points {

namespace internal {

/**
 *   @brief Interface to access X coordinate value
 */
template <typename real>
requires std::is_floating_point_v<real>
struct RequireXaccessInterface {
  /**
   *  @brief Return value of X coordinate
   *
   *  @return X coordinate for reading
   */
  [[nodiscard]] inline constexpr real x() const noexcept;

  /**
   *  @brief Return reference to X coordinate for modifying purposes
   *
   *  @return reference to X coordinate for writing
   */
  [[nodiscard]] inline constexpr real &x() noexcept;
};

/**
 *   @brief Interface to access Y coordinate value
 */
template <typename real>
requires std::is_floating_point_v<real>
struct RequireYaccessInterface {
  /**
   *  @brief Return value of Y coordinate
   *
   *  @return Y coordinate for reading
   */
  [[nodiscard]] inline constexpr real y() const noexcept;

  /**
   *  @brief Return reference to Y coordinate for modifying purposes
   *
   *  @return reference to Y coordinate for writing
   */
  [[nodiscard]] inline constexpr real &y() noexcept;
};

/**
 *   @brief Interface to access Z coordinate value
 */
template <typename real>
requires std::is_floating_point_v<real>
struct RequireZaccessInterface {
  /**
   *  @brief Return value of Z coordinate
   *
   *  @return Z coordinate for reading
   */
  [[nodiscard]] inline constexpr real z() const noexcept;

  /**
   *  @brief Return reference to Z coordinate for modifying purposes
   *
   *  @return reference to Z coordinate for writing
   */
  [[nodiscard]] inline constexpr real &z() noexcept;
};

template <class T> struct NotRequired {};

template <bool b, class T> using AddInterfaceIf = typename std::conditional<b, T, internal::NotRequired<T>>::type;
} // namespace internal

/**
 *  @brief A base class for point implementations
 * \See internal::RequireXaccessInterface
 * \See internal::RequireYaccessInterface
 * \See internal::RequireZaccessInterface
 */
template <std::size_t DIM, typename type = float>
requires std::is_floating_point_v<type> &&((DIM == 2u) || (DIM == 3u)) /**/
    class PointD : public internal::AddInterfaceIf<DIM >= 1u, internal::RequireXaccessInterface<type>>,
                   internal::AddInterfaceIf<DIM >= 2u, internal::RequireYaccessInterface<type>>,
                   internal::AddInterfaceIf<DIM >= 3u, internal::RequireZaccessInterface<type>> {
protected:
  static constexpr std::size_t const dim = DIM;
  /// Define the dimension
  enum class Dimension : std::size_t;

private:
  using Point = PointD<DIM, type>;

public:
  using real = type;

  /**
   *  @brief The default constructor. Coordinated are zeroed
   */
  inline constexpr PointD() noexcept = default;

  /**
   *  @brief The destructor
   */
  inline constexpr ~PointD() noexcept = default;

public:
  /**
   *  @brief Multiply coordinate values of point with 'scale'
   *
   *  @param scale multiplication factor
   *  @return new point that is multiplied
   */
  inline constexpr PointD operator*(real const scale) const noexcept;

  /**
   *  @brief Dot product of two points
   *
   *  @param p another point
   *  @return Dot product of two points
   */
  inline constexpr real operator*(PointD const p) const noexcept;

  /**
   *  @brief Divide coodinate vales of point by 'div'
   *
   *  @param div divider
   *  @return the divided point as std::optional
   *  @return std::nullopt, if divided by zero or near zero
   */
  [[nodiscard]] inline constexpr std::optional<PointD> operator/(real const div) const noexcept;

  /**
   *  @brief Calculate substraction between this and another point 'p'
   *
   *  @param p point to be subracted from 'this'
   *  @return a new point containing the substraction
   */
  [[nodiscard]] inline constexpr PointD operator-(PointD const p) const noexcept;

  /**
   *  @brief Calculate addition between this and point 'p'
   *
   *  @param p point to be added to 'this'
   *  @return a new point containing the addition
   */
  [[nodiscard]] inline constexpr PointD operator+(PointD const p) const noexcept;

  /**
   *  @brief Add point 'p' to this point
   *
   *  @param p point to be added to 'this'
   *  @return 'this' point
   */
  inline constexpr PointD &operator+=(PointD const p) noexcept;

  /**
   *  @brief Calculate trace of point or sum of coordinate values. This is different than 1-norm:
   * \f$\stackrel[i=1]{n}{\sum}\left|x_{i}\right|\f$
   *
   *  @return sum of coordinates
   */
  [[nodiscard]] inline constexpr real trace() const noexcept;

  /**
   *  @brief Calculate squared length of point as vector
   *
   *  @return squared length
   */
  [[nodiscard]] inline constexpr real lengthSquared() const noexcept;

  /**
   *  @brief Calcualate length of point as vector
   *
   *  @return length
   */
  [[nodiscard]] inline constexpr real length() const noexcept;

  /**
   *  @brief Calculate distance between this point and point 'p'
   *
   *  @param p the point calcualte distance to
   *  @return distance between points
   */
  [[nodiscard]] inline constexpr real distance(PointD const p) const noexcept;

  /**
   *  @brief Normalize lenth of point (as a vector) to unit length
   *
   *  @return a point, which length is one as std::optional
   *  @return std::nullopt, if length of point is zero or near zero
   */
  [[nodiscard]] inline constexpr std::optional<PointD> normalize() const noexcept;

  /**
   *  @brief Return order or number of coordinate values point have
   *
   *  @return DIM
   */
  [[nodiscard]] static inline consteval std::size_t getOrder() noexcept;

protected:
  /// The coordinate values
  std::array<real, DIM> coords_;
};

} // namespace curve::points

#endif
