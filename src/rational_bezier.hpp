//
// (C) Matti Lehtonen 2023
//

/**
 *  @file rational_bezier.hpp This file contains implementation of rational bezier class.
 */

#ifndef RATIONAL_BEZIER_H
#define RATIONAL_BEZIER_H

#include "validate.hpp"

namespace curve {
namespace bezier {
namespace rational {

/**
 *  @brief The class for rational bezier curves with management of control points of type 'CP', either Point2<> or
 * Point3<>. Delegates calls to @ref curve::bezier::rational::ValidateRational
 */
template <class CP>
requires std::is_same_v<CP, ControlPoint<typename CP::Point>>
class Rational : public DelegationInterface<CP> {
private:
  using Interface = DelegationInterface<CP>;

public:
  using ControlPoint = typename Interface::ControlPoint;
  using ConstControlPoint = ControlPoint const;
  using Point = typename Interface::Point;
  using ConstPoint = Point const;
  using real = typename Interface::Point::real;

  using ControlPointContainer = typename std::vector<ControlPoint>;
  using ControlPointSpan = typename std::span<ConstControlPoint>;
  using Validator = ValidateRational<ControlPoint>;

  using point_or_issue = typename Interface::point_or_issue;
  using real_or_issue = typename Interface::real_or_issue;
  using vector_of_points_or_issue = typename Interface::vector_of_points_or_issue;

  /**
   *  @brief the default constructor
   */
  inline constexpr Rational() = default;

  /**
   *  @brief Constructor with vector of control points copied to the curve.
   *  NOTE: cumulative weight of control points must NOT be near zero.
   *
   *  @param cp The control points of 'CP to be copied to this curve
   */
  inline constexpr Rational(ControlPointContainer const &cp) {
    this->controlPoints_.reserve(cp.size()); // Might throw
    this->controlPoints_.insert(this->controlPoints_.begin(), cp.cbegin(), cp.cend());
  }

  /**
   *  @brief Constructor with vector of control points moved to the curve.
   *  NOTE: cumulative weight of control points must not be near zero.
   *
   *  @param cp The control points of 'CP to be moved to this curve
   */
  inline constexpr Rational(ControlPointContainer &&cp) noexcept { this->controlPoints_ = std::move(cp); }

  /**
   *  @brief Destructor for rational bezier curve
   */
  inline constexpr ~Rational() noexcept { this->controlPoints_ = ControlPointContainer(); }

  /**
   *  @brief Output bezier curve 'data' to stream 'out'
   *
   *  @param out the output stream
   *  @param data The rational bezier curve to be outputted
   *  @return the 'out' stream
   */
  inline friend std::ostream &operator<<(std::ostream &out, Rational const &data) {
    out << "{";
    out << data.controlPoints_;
    out << "}";
    return out;
  }

public:
  /**
   *  @brief Calculate the point from the curve at point 'u'
   *
   *  @param u a position from curve, range [0;1]
   *  @return a point from curve, depending type of 'CP', either 2D or 3D, as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr point_or_issue pointAt(real const u) const noexcept {
    Validator const validator(this->controlPointContainerAsSpan());
    return validator.pointAt(u);
  }

  /**
   *  @brief Calculate the length of the curve
   *
   *  @return the length of curve as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr real_or_issue curveLength() const noexcept {
    Validator const validator(this->controlPointContainerAsSpan());
    return validator.curveLength();
  }

  /**
   *  @brief Generate an approximation of curve as a line string
   *
   *  @return a vector containing vertices between line segments as std::variant
   *  @return a vector containning a single vertex, if curve has only one control point as std::optional
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr vector_of_points_or_issue asLinestring() const noexcept {
    Validator const validator(this->controlPointContainerAsSpan());
    return validator.asLineString();
  }

  /**
   *  @brief Calculate velocity from curve at point 'u'
   *
   *  @param u a position from curve, range [0;1]
   *  @return velocity at point from curve, depending type of 'CP', either 2D or 3D, as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr point_or_issue velocityAt(real const u) const noexcept {
    Validator const validator(this->controlPointContainerAsSpan());
    return validator.velocityAt(u);
  }

  /**
   *  @brief Calculate speed from curve at point 'u'
   *
   *  @param u a position from curve, range [0;1]
   *  @return speed at point from curve, as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr real_or_issue speedAt(real const u) const noexcept {
    Validator const validator(this->controlPointContainerAsSpan());
    return validator.speedAt(u);
  }

  /**
   *  @brief Calculate unit length tangent from curve at point 'u'
   *
   *  @param u a position from curve, range [0;1]
   *  @return tangent (unit length) at point from curve, depending type of 'CP', either 2D or 3D, as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr point_or_issue tangentAt(real const u) const noexcept {
    Validator const validator(this->controlPointContainerAsSpan());
    return validator.tangentAt(u);
  }

  /**
   *  This problem is actually quite challenging than usual find global minimum of function. Normally the minimum of
   * function D_p(u) is located at either: a) D_p'(u) = 0 b) ends of search range, u = {0,1} c) discontinuation
   * locations of D_p(u), if such exists
   *
   *  How ever, now D_p(u) is expected to have multiple local minimums (think curve like snake) and in an extreme
   * case D_p(u) can be constant (think curve like arc of circle, where p is origo of that circle).
   *  =-> D_p'(u) = D_p"(u) = 0 or no global minimum exists as a single point at curve.
   *
   *  See e.g. of research related to the topic
   *  Xiao-Diao Chen, Yin Zhou, Zhenyu Shu, Hua Su, Jean-Claude Paul. Improved Algebraic Algorithm On Point
   * Projection For Bézier Curves. Proceedings of the Second International Multi-Symposiums on Computer and
   * Computational Sciences (IMSCCS 2007), The University of Iowa, Iowa City, Iowa, USA, Aug 2007, Iowa, United
   * States. pp.158-163, 10.1109/IMSCCS.2007.17. inria-00518379
   *
   *  https://pomax.github.io/bezierinfo/#projections
   *  https://pomax.github.io/bezierinfo/#extremities
   *
   *  NOTE!:
   *  The current implementation use Newton-Raphson method to find locations, where derivate of distance (squared)
   * between given point and curve is zero or df(u)/du = 0.
   *
   *  @brief Find any point (potentially of many), which is closest to given point 'p'. First generates a rought
   * approximation of curve as line string with 2N-1 vertices. Then use these 2N-1 vertices as initial guesses for
   * Newton-Raphson method to find locations, where either the function reaches zero f(x) = 0 or the derivate f'(x)
   * goes to (near) zero. Finally selects the point with the smallest distance to 'p', if any. For function f(x) is
   * used a derivate of distance squared between curve and point 'p'. For function f'(x) is used an approximation of
   * 2nd derivative of distace squared between curve and point 'p'. See e.g.
   * https://math.hmc.edu/calculus/hmc-mathematics-calculus-online-tutorials/single-variable-calculus/limit-definition-of-the-derivative/
   *
   *  @param p The point near (or on) curve
   *  @return one of points (if many) that is the closest one to the point 'p' as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline point_or_issue closestCurvePointFor(ConstPoint p) const noexcept {
    Validator const validator(this->controlPointContainerAsSpan());
    return validator.closestCurvePointFor(p);
  }

  /**
   *  @brief Return number of control point defined with the curve
   *
   *  @return Number of control point defined for the curve
   */
  [[nodiscard]] inline constexpr std::size_t numberOfControlPoints() const noexcept {
    return this->controlPoints_.size();
  }

  /**
   *  @brief Get a reference to control point at index 'index' for modifying purposes.
   * NOTE1: this reference may get invalidated, if allocation of container of control points is changed.
   * NOTE2: If weight is modified so that sum of all weights goes near zero, the curve shall be invalid
   *
   *  @param index the index of control point, valid range [ 0 ; \<number of control points\> [
   *  @return Reference to control point to modify as std::optional
   *  @return std::nullopt, if input is invalid
   */
  [[nodiscard]] inline constexpr std::optional<std::reference_wrapper<ControlPoint>> const
  getControlPoint(std::size_t const index) noexcept {
    if (index < this->numberOfControlPoints()) {
      return std::make_optional(std::ref(this->controlPoints_.at(index)));
    } else {
      return std::nullopt;
    }
  }

  /**
   *  @brief Get a const reference to control point at index 'index' for reading purposes
   * NOTE: this reference may get invalidated, if allocation of container of control points is changed.
   *
   *  @param index the index of control point, valid range [ 0 ; \<number of control points> \[
   *  @return Const control point for read access as std::optional
   *  @return std::nullopt, if input is invalid
   */
  [[nodiscard]] inline constexpr std::optional<ControlPoint const> const
  getControlPoint(std::size_t const index) const noexcept {
    if (index < this->numberOfControlPoints()) {
      return std::make_optional(std::cref(this->controlPoints_.at(index)));
    } else {
      return std::nullopt;
    }
  }

  /**
   *  @brief Add a new control point 'cp' to the curve at index 'index'
   *
   *  @param index the index of control point, valid range [ 0 ; \<number of control points\> [
   *  @param cp the added control point of type 'CP'
   *  @return true, if the index is within range [ 0 ; \<number of control points\> [
   *  @return false otherwise
   */
  [[nodiscard]] inline constexpr bool addControlPoint(std::size_t const index, ControlPoint const cp) {
    if (index <= this->numberOfControlPoints()) {
      try {
        (void)this->controlPoints_.emplace(this->controlPoints_.begin() + index, cp); // Might throw
      }

      catch (std::exception const &) {
        return false;
      }

      return true;
    } else {
      return false;
    }
  }

  /**
   *  @brief Remove a control point from curve at index 'index'
   *
   *  @param index the index of control point, valid range [ 0 ; \<number of control points\> [
   *  @return true, if the index is within range [ 0 ; \<number of control points\> [
   *  @return false otherwise
   */
  [[nodiscard]] inline constexpr bool removeControlPoint(std::size_t const index) noexcept {
    if (index >= this->controlPoints_.size()) {
      return false;
    }
    (void)this->controlPoints_.erase(this->controlPoints_.begin() + index);
    return true;
  }

  /**
   *  @brief Remove all control points from curve, allocation of container for control points is not changed.
   */
  inline constexpr void removeAllControlPoints() noexcept { this->controlPoints_.clear(); }

  /**
   *  @brief Return the current capacity of container for control points
   *
   *  @return Number of control points that can be stored to the container without increasing allocation for
   * container
   */
  inline constexpr std::size_t capacity() const noexcept { return this->controlPoints_.capacity(); }

  /**
   *  @brief Reserve 'n' number of control points at container. Allocation is changed, if 'n' is greater than the
   * current capacity
   *
   *  @param n reserve space for 'n' control points
   */
  inline constexpr bool reserve(std::size_t const n) {
    try {
      this->controlPoints_.reserve(n); // Might throw
    }

    catch (std::exception const &) {
      return false;
    }

    return true;
  }

  /**
   *  @brief Match capacity of container of control points with the number of control points at container
   */
  inline constexpr void slim() noexcept { this->controlPoints_.shrink_to_fit(); }

protected:
  [[nodiscard]] inline constexpr ControlPointSpan controlPointContainerAsSpan() const noexcept {
    return ControlPointSpan(this->controlPoints_.cbegin(), this->controlPoints_.size());
  }

protected:
  ControlPointContainer controlPoints_; //< The container for control points
};

} // namespace rational
} // namespace bezier
} // namespace curve

#endif
