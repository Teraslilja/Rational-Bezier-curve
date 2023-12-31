//
// (C) Matti Lehtonen 2023
//

/**
 *  @file validation.hpp This file contains module interface for validation of rational bezier
 */

#ifndef VALIDATION_H
#define VALIDATION_H

#include "calculation.hpp"
#include "validation_interface.hpp"

namespace curve::bezier::rational {

/**
 *  @brief The class to validate rational bezier curves with control point of type 'CP', using either Point2<> or
 * Point3<> and if curve is valid, calculate points from it.
 */
template <class CP>
requires std::is_same_v<CP, ControlPoint<typename CP::Point>>
class ValidateRational : private curve::bezier::rational::ValidationInterface<CP> {
private:
  using Interface = curve::bezier::rational::ValidationInterface<CP>;
  using Calculator = curve::bezier::rational::CalculateRational<CP>;

public:
  using ControlPoint = CP;
  using Point = typename ControlPoint::Point;
  using real = typename Point::real;

  using ControlPointSpan = typename std::span<ControlPoint const>;

  using point_or_issue = typename Interface::point_or_issue;
  using real_or_issue = typename Interface::real_or_issue;
  using vector_of_points_or_issue = typename Interface::vector_of_points_or_issue;

  /**
   *  @brief Move constructor for object to validate and calculate rational bezier curve from given control points
   *
   *  @param controlPoints read-only span from actual container of control points
   */
  inline constexpr ValidateRational(ControlPointSpan const &&controlPoints) noexcept
      : calculator(std::move(controlPoints)) {}

  /**
   *  @brief A default destructor
   */
  inline constexpr ~ValidateRational() noexcept = default;

private:
  /**
   *  @brief Check, if curve parameter 'u' is in valid range [0;1]
   *
   *  @param u curve parameter, valid range is [0;1]
   *  @return std::nullopt, if parameter 'u' is valid
   *  @return ValidityIssue::ISSUE_U_IS_INVALID as std::optional otherwise
   */
  [[nodiscard]] static inline constexpr std::optional<ValidityIssue> isValidU(real const u) noexcept {
    return ((u >= Calculator::U_MIN) && (u <= Calculator::U_MAX))
               ? std::nullopt
               : std::make_optional(ValidityIssue::ISSUE_U_IS_INVALID);
  }

  /**
   *  @brief Check, if curve has enough (>= 'count') control points
   *
   *  @param count Wanted minimum amount of control points
   *  @return std::nullopt, if enough control points
   *  @return ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS as std::optional otherwise
   */
  [[nodiscard]] inline constexpr std::optional<ValidityIssue>
  hasEnoughControlPoints(std::size_t const count) const noexcept {
    return (this->calculator.numberOfControlPoints() >= count)
               ? std::nullopt
               : std::make_optional(ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS);
  }

  /**
   *   @brief Is weight of control point valid
   *
   *   @param cp control point
   *  @return true, if weight is valid
   *  @return false otherwise
   */
  [[nodiscard]] static inline constexpr bool isWeightOfControlPointValid(ControlPoint const &cp) noexcept {
    return std::isfinite(cp.w());
  }

  /**
   *   @brief Are weights of control points valid
   *
   *  @return std::nullopt, if weights are valid
   *  @return ValidityIssue::ISSUE_BAD_CONTROLPOINT_WEIGHT as std::optional otherwise
   */
  [[nodiscard]] inline constexpr std::optional<ValidityIssue> hasWeightsOfControlPointsValid() const noexcept {
    ControlPointSpan const &cps = this->calculator.getSpan();
    for (auto const &cp : cps) {
      bool const valid = isWeightOfControlPointValid(cp);
      if (!valid) {
        return ValidityIssue::ISSUE_BAD_CONTROLPOINT_WEIGHT;
      }
    }
    return std::nullopt;
  }

  /**
   *  @brief Is coodenate value of point valid
   *
   *  @return true, if coordinate value is valid
   *  @return false otherwise
   */
  [[nodiscard]] static inline constexpr bool isValidCoordinateValueOfPoint(real const v) noexcept {
    return std::isfinite(v);
  }

  /**
   *  @brief Are points of control points valid
   *
   *  @param cp control point
   *  @return true, if point is valid
   *  @return false otherwise
   */
  [[nodiscard]] static inline constexpr bool isValidPointOfControlPoint(ControlPoint const &cp) noexcept {
    Point const &p = cp.p();
    bool isValid = (p.getOrder() > 0u);
    if constexpr (p.getOrder() >= 1u) {
      isValid = isValid && isValidCoordinateValueOfPoint(p.x());
    } else if constexpr (p.getOrder() >= 2u) {
      isValid = isValid && isValidCoordinateValueOfPoint(p.y());
    } else if constexpr (p.getOrder() >= 3u) {
      isValid = isValid && isValidCoordinateValueOfPoint(p.z());
    }

    return isValid;
  }

  /**
   *  @brief Are points of control points valid
   *
   *  @return std::nullopt, if points are valid
   *  @return ValidityIssue::ISSUE_BAD_POINT as std::optional otherwise
   */
  [[nodiscard]] inline constexpr std::optional<ValidityIssue> hasValidPoints() const noexcept {
    ControlPointSpan const &cps = this->calculator.getSpan();
    for (auto const &cp : cps) {
      bool const isValid = isValidPointOfControlPoint(cp);
      if (!isValid) {
        return ValidityIssue::ISSUE_BAD_POINT;
      }
    }
    return std::nullopt;
  }

  /**
   *  @brief Check, if sum of curve's weights of control points isn't too close to zero
   *
   *  @param u curve parameter, valid range is [0;1]
   *  @return std::nullopt, if sum of weights differs from zero more than real::epsilon
   *  @return ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS as std::optional otherwise
   */
  [[nodiscard]] inline constexpr std::optional<ValidityIssue> hasValidSumOfWeights(real const u) const noexcept {
    real const sum_w = [u, this]() -> real {
      ControlPointSpan const &cps = this->calculator.getSpan();
      std::size_t const n = cps.size() - 1u;
      real sum_w = real(0);
      for (std::size_t i = 0u; auto const &cp : cps) {
        real const w = curve::bezier::utilities::BernsteinPolynomials<real>::B(i, n, u) * cp.w();
        sum_w += w;

        ++i;
      }
      return sum_w;
    }();

    return (std::abs(sum_w) >= std::numeric_limits<real>::epsilon())
               ? std::nullopt
               : std::make_optional(ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS);
  }

  /**
   *  @brief Check, if curve's control point configuration is valid (enough control points and okay weights)
   *
   *  @param count Wanted minimum amount of control points
   *  @return std::nullopt, if enough control points
   *  @return ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS as std::optional otherwise
   */
  [[nodiscard]] inline constexpr std::optional<ValidityIssue> isValid(std::size_t const count) const noexcept {
    std::optional<ValidityIssue> issue;

    issue = this->hasEnoughControlPoints(count);
    if (issue.has_value()) {
      return issue.value();
    }

    issue = hasValidPoints();
    if (issue.has_value()) {
      return issue.value();
    }

    issue = hasWeightsOfControlPointsValid();
    if (issue.has_value()) {
      return issue.value();
    }

    return this->hasValidSumOfWeights(real(0.5));
  }

  /**
   *  @brief Check, if curve's control point configuration is valid (enough control points and okay weights) and
   * selecting a valid point (parameter u) from curve
   *
   *  @param count Wanted minimum amount of control points
   *  @param u curve parameter, valid range is [0;1]
   *  @return std::nullopt, if sum of weights differs from zero more than real::epsilon
   *  @return ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS as std::optional otherwise
   */
  [[nodiscard]] inline constexpr std::optional<ValidityIssue> isValid(std::size_t const count,
                                                                      real const u) const noexcept {
    std::optional<ValidityIssue> issue;
    issue = this->hasEnoughControlPoints(count);
    if (issue.has_value()) {
      return issue.value();
    }

    issue = this->isValidU(u);
    if (issue.has_value()) {
      return issue.value();
    }

    issue = hasValidPoints();
    if (issue.has_value()) {
      return issue.value();
    }

    issue = hasWeightsOfControlPointsValid();
    if (issue.has_value()) {
      return issue.value();
    }

    return this->hasValidSumOfWeights(u);
  }

  /**
   *  @brief Check, if curve's control point configuration is valid (enough control points and okay weights) and
   * selecting a valid point (parameter u) from curve
   *
   *  @param count Wanted minimum amount of control points
   *  @param us curve parameters, valid range is [0;1]
   *  @return std::nullopt, if sum of weights differs from zero more than real::epsilon
   *  @return ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS as std::optional otherwise
   */
  [[nodiscard]] inline constexpr std::optional<ValidityIssue> isValid(std::size_t const count,
                                                                      std::span<real> const us) const noexcept {
    std::optional<ValidityIssue> issue;
    issue = this->hasEnoughControlPoints(count);
    if (issue.has_value()) {
      return issue.value();
    }

    issue = hasValidPoints();
    if (issue.has_value()) {
      return issue.value();
    }

    issue = hasWeightsOfControlPointsValid();
    if (issue.has_value()) {
      return issue.value();
    }

    for (real const u : us) {
      issue = this->isValidU(u);
      if (issue.has_value()) {
        return issue.value();
      }

      issue = this->hasValidSumOfWeights(u);
      if (issue.has_value()) {
        return issue.value();
      }
    }

    return std::nullopt;
  }

public:
  /**
   *  @brief Calculate the point from the curve at point 'u'
   *
   *  @param u a position from curve, range [0;1]
   *  @return a point from curve, depending type of 'CP', either 2D or 3D, as std::variant
   *  @return Issue ID as std::variant, if curve is invalid
   */
  [[nodiscard]] inline constexpr point_or_issue pointAt(real const u) const noexcept {
    std::optional<ValidityIssue> const issue = this->isValid(1u, u);
    if (issue.has_value()) {
      return issue.value();
    }

    return this->calculator.C(u);
  }

  /**
   *  @brief Calculate multible points from curve
   *
   *  @param us a span of positions from curve, all within range [0;1]
   *  @return vector of newly calculated vertives (as std::variant> or
   *  @return Issue number (as std::variant), if curve is invalid or problems with memory (re)allocation
   */
  [[nodiscard]] constexpr vector_of_points_or_issue pointsAt(std::span<real> const us) const {
    std::optional<ValidityIssue> const issue = this->isValid(1u, us);
    if (issue.has_value()) {
      return issue.value();
    }

    std::vector<Point> points;
    try {
      points.reserve(us.size()); // Might throw
    }

    catch (std::exception const &) {
      return ValidityIssue::ISSUE_OUT_OF_HEAP_MEMORY;
    }

    for (real const u : us) {
      points.emplace_back(this->C(u));
    }

    return std::move(points);
  }

  /**
   *  @brief Calculate the length of the curve
   *
   *  @return the length of curve as std::variant
   *  @return Issue numner, if curve is invalid
   */
  [[nodiscard]] inline constexpr real_or_issue curveLength() const noexcept {
    std::optional<ValidityIssue> const issue = this->isValid(1u);
    if (issue.has_value()) {
      return issue.value();
    }

    if (this->calculator.numberOfControlPoints() < 2u) {
      // point as a degenerate curve
      return real(0);
    } else {
      return this->calculator.curveLength();
    }
  }

  /**
   *  @brief Generate an approximation of curve as a line string
   *
   *  @return a vector containing vertices between line segments as std::variant
   *  @return a vector containning a single vertex, if curve has only one control point as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr vector_of_points_or_issue asLineString() const noexcept {
    std::optional<ValidityIssue> const issue = this->isValid(1u);
    if (issue.has_value()) {
      return issue.value();
    }

    if (this->calculator.numberOfControlPoints() < 2u) {
      // point as a degenerate curve
      std::vector<Point> single_point(1u, this->calculator.getSpan().front().p());
      return single_point;
    } else {
      std::vector<Point> lineSegments = this->calculator.asLinestring();
      return lineSegments;
    }
  }

  /**
   *  @brief Calculate velocity from curve at point 'u'
   *
   *  @param u a position from curve, range [0;1]
   *  @return velocity at point from curve, depending type of 'CP', either 2D or 3D, as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr point_or_issue velocityAt(real const u) const noexcept {
    std::optional<ValidityIssue> const issue = this->isValid(2u, u);
    if (issue.has_value()) {
      return issue.value();
    }

    return this->calculator.dC(u);
  }

  /**
   *  @brief Calculate speed from curve at point 'u'
   *
   *  @param u a position from curve, range [0;1]
   *  @return speed at point from curve, as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr real_or_issue speedAt(real const u) const noexcept {
    std::optional<ValidityIssue> const issue = this->isValid(2u, u);
    if (issue.has_value()) {
      return issue.value();
    }

    Point const velocity = this->calculator.dC(u);
    return velocity.length();
  }

  /**
   *  @brief Calculate unit length tangent from curve at point 'u'
   *
   *  @param u a position from curve, range [0;1]
   *  @return tangent (unit length) at point from curve, depending type of 'CP', either 2D or 3D, as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   *  @return vector of {0,0,0}, if cannot normalize tangent (shouldn't be possible)
   */
  [[nodiscard]] inline constexpr point_or_issue tangentAt(real const u) const noexcept {
    std::optional<ValidityIssue> const issue = this->isValid(2u, u);
    if (issue.has_value()) {
      return issue.value();
    }

    Point const velocity = this->calculator.dC(u);

    // It should be impossible NOT to normalize tangent ..
    std::optional<Point const> const has_unit_length = velocity.normalize();
    // .. but, if that happens, return vector with zero length
    constexpr Point const zero_length = Point();
    return has_unit_length.value_or(zero_length);
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
  [[nodiscard]] inline constexpr point_or_issue closestCurvePointFor(Point const p) const noexcept {
    std::optional<ValidityIssue> const issue = this->isValid(1u);
    if (issue.has_value()) {
      return issue.value();
    }

    if (this->calculator.numberOfControlPoints() == 1u) {
      // point as a degenerate curve
      return this->calculator.getSpan()[0u].p();
    } else {
      return this->calculator.findNearestPointFor(p);
    }
  }

private:
  Calculator const calculator;
};

} // namespace curve::bezier::rational

#endif
