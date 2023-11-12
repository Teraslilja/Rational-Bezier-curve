//
// (C) Matti Lehtonen 2023
//

/**
 *  @file rational_bezier.hpp This file contains implementation of rational bezier class.
 */

#ifndef RATIONAL_BEZIER_H
#define RATIONAL_BEZIER_H

#include "bernstein_polynomials.hpp"
#include "control_point.hpp"
#include "newton_raphson.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <optional>
#include <span>
#include <type_traits>
#include <variant>
#include <vector>

namespace curve {
namespace bezier {
namespace rational {
namespace internal {

/**
 *  @brief The class to calculate rational bezier curves with control point of type 'CP', using either Point2<> or
 * Point3<>
 */
template <class CP>
requires std::is_same_v<CP, const ControlPoint<typename CP::point>>
class CalculateRational {
public:
  using ControlPoint = CP; //< This class uses only read-only control points
  using point = typename ControlPoint::point;
  using real = typename point::real;
  using span = typename std::span<ControlPoint const>;

  static real constexpr U_MIN = real(0); //< Minimum valid value of u
  static real constexpr U_MAX = real(1); //< Maximum valid value of u

  /**
   *  @brief Move constructor for object to calculate rational bezier curve from given control points
   *
   *  @param controlPoints read-only span from actual container of control points
   */
  inline constexpr CalculateRational(span const &&controlPoints) noexcept
      : controlPoints_(std::move(controlPoints)) {}

  /**
   *  @brief A default destructor
   */
  inline constexpr ~CalculateRational() = default;

  /**
   *  @brief Function C(u) to calculate a point at 'u'. @see pointAt()
   *
   *  @param u a position from curve, range [0;1]
   *  @return a point from curve, depending type of 'CP', either 2D or 3D
   */
  [[nodiscard]] constexpr point C(real const u) const noexcept {
    /**
       \f$C(u)=\frac{f(u)}{g(u)}\f$, where
       \f$f(u)=\underset{i=0}{\overset{n}{\sum}}B_{i,n}(u)w_{i}\mathbf{P}_{i}\f$
       and \f$g(u)=\underset{i=0}{\overset{n}{\sum}}B_{i,n}(u)w_{i}\f$
    */

    // Calculate f(u) and g(u)
    auto const [f_u, g_u] = [u](span const cp) -> std::tuple<point, real> {
      std::size_t const n = cp.size() - 1u;
      point sum_wp;
      real sum_w = real(0);
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::B(i, n, u) * cp[i].w();
        sum_wp += (cp[i].p() * w);
        sum_w += w;
      }
      return {sum_wp, sum_w};
    }(this->controlPoints_);

    std::optional<point> const C_u = f_u / g_u;
    return C_u.value();
  }

  /**
   *  @brief Return number of control points
   *
   *  @return Number of control points
   */
  [[nodiscard]] inline constexpr std::size_t numberOfControlPoints() const noexcept {
    return this->controlPoints_.size();
  }

  /**
   *  @brief Return the span of control points
   *
   *  @return reference to the span of control points
   */
  [[nodiscard]] inline constexpr span const &getSpan() const noexcept { return this->controlPoints_; }

  /**
   *  @brief Function C'(u) to calculate a derivate at 'u'. @see velocityAt()
   *
   *  @param u a position from curve, range [0;1]
   *  @return a point from curve, depending type of 'CP', either 2D or 3D
   */
  [[nodiscard]] constexpr point dC(real const u) const noexcept {
    /**
       \f$C(u)=\frac{f(u)}{g(u)}\f$, where
       \f$f(u)=\underset{i=0}{\overset{n}{\sum}}B_{i,n}(u)w_{i}\mathbf{P}_{i}\f$
       and \f$g(u)=\underset{i=0}{\overset{n}{\sum}}B_{i,n}(u)w_{i}\f$

       \f$\Rightarrow C'(u)=\frac{f'(u)g(u)-f(u)g'(u)}{g(u)^{2}}\f$, where
       \f$f'(u)=\underset{i=0}{\overset{n}{\sum}}B'_{i,n}(u)w_{i}\mathbf{P}_{i}\f$
       and \f$g'(u)=\underset{i=0}{\overset{n}{\sum}}B'_{i,n}(u)w_{i}\f$,
       where \f$B'_{i,n}(u)=n\left(B_{i-1,n-1}(u)-B_{i,n-1}(u)\right)\f$
    */

    // Calculate f(u) and g(u)
    auto const [f_u, g_u] = [u](span const &cp) -> std::tuple<point, real> {
      std::size_t const n = cp.size() - 1u;
      real sum_w = real(0);
      point sum_wp;
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::B(i, n, u) * cp[i].w();
        sum_w += w;
        sum_wp += (cp[i].p() * w);
      }
      return {sum_wp, sum_w};
    }(this->controlPoints_);

    // Calculate f'(u) and g'(u)
    auto const [df_u, dg_u] = [u](span const &cp) -> std::tuple<point, real> {
      std::size_t const n = cp.size() - 1u;
      point sum_wp;
      real sum_w = real(0);
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::dB(i, n, u) * cp[i].w();
        sum_wp += cp[i].p() * w;
        sum_w += w;
      }
      return {sum_wp, sum_w};
    }(this->controlPoints_);

    // C'(u) = (f'(u) * g(u) - f(u) * g'(u)) / g(u)^2
    std::optional<point> const dC_u = ((df_u * g_u) - (f_u * dg_u)) / (g_u * g_u);
    return dC_u.value();
  }

  /**
   *  @brief Function C"(u) to calculate a second derivate at 'u'. @see
   * findNearestPointFor()
   *
   *  @param controlPoints the control points for rational belizer
   *  @param u a position from curve, range [0;1]
   *  @return a point from curve, depending type of 'CP', either 2D or 3D
   */
  [[nodiscard]] constexpr point d2C(span const controlPoints, real const u) noexcept {
    /**
       \f$C(u)=\frac{f(u)}{g(u)}\f$, where
       \f$f(u)=\underset{i=0}{\overset{n}{\sum}}B_{i,n}(u)w_{i}\mathbf{P}_{i}\f$
       and \f$g(u)=\underset{i=0}{\overset{n}{\sum}}B_{i,n}(u)w_{i}\f$

       \f$\Rightarrow
       C"(u)=\frac{g^{2}(u)f"(u)-g(u)\left(2f'(u)g'(u)+f(x)g"(x)\right)+2f(x)g'^{2}(x)}{g^{3}(u)}\f$
       \f$f'(u)=\underset{i=0}{\overset{n}{\sum}}B'_{i,n}(u)w_{i}\mathbf{P}_{i}\f$,
       \f$g'(u)=\underset{i=0}{\overset{n}{\sum}}B'_{i,n}(u)w_{i}\f$,
       \f$f"(u)=\underset{i=0}{\overset{n}{\sum}}B"_{i,n}(u)w_{i}\mathbf{P}_{i}\f$,
       and \f$g"(u)=\underset{i=0}{\overset{n}{\sum}}B"_{i,n}(u)w_{i}\f$.
    */

    // Calculate f(u) and g(u)
    auto const [f_u, g_u] = [u](span const &cp) -> std::tuple<point, real> {
      std::size_t const n = cp.size() - 1u;
      real sum_w = real(0);
      point sum_wp;
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::B(i, n, u) * cp[i].w();
        sum_w += w;
        sum_wp += (cp[i].p() * w);
      }
      return {sum_wp, sum_w};
    }(controlPoints);

    // Calculate f'(u) and g'(u)
    auto const [df_u, dg_u] = [u](span const &cp) -> std::tuple<point, real> {
      std::size_t const n = cp.size() - 1u;
      point sum_wp;
      real sum_w = real(0);
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::dB(i, n, u) * cp[i].w();
        sum_wp += cp[i].p() * w;
        sum_w += w;
      }
      return {sum_wp, sum_w};
    }(controlPoints);

    // Calculate f"(u) and g"(u)
    auto const [d2f_u, d2g_u] = [u](span const &cp) -> std::tuple<point, real> {
      std::size_t const n = cp.size() - 1u;
      point sum_wp;
      real sum_w = real(0);
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::d2B(i, n, u) * cp[i].w();
        sum_wp += cp[i].p() * w;
        sum_w += w;
      }
      return {sum_wp, sum_w};
    }(controlPoints);

    // C"(u) = ( g'(u)^2 * f"(u) - g(u) * (2 * f'(u)*g'(u) + f(u) * g"(u)) + 2 *
    // f(u) * g'(u)^2) / g(u)^3
    std::optional<point> const d2C_u =
        (d2f_u * (dg_u * dg_u) - (df_u * (real(2) * dg_u) + f_u * d2g_u) * g_u + f_u * (real(2) * dg_u * dg_u)) /
        (g_u * g_u * g_u);
    return d2C_u.value();
  }

  /**
   *  @brief Calculate the length of the curve
   *
   *  @return the length of curve
   */
  [[nodiscard]] constexpr real length() const noexcept {
    real constexpr REQUIRED_SEGMENT_LENGTH_ACCURACY = real(1e-6);

    real length = real(0);
    auto const cumulate_segment_lengths = [&length](real const segment_length,
                                                    std::pair<real, point> const segment_end) -> void {
      (void)segment_end;
      length += segment_length;
    };

    this->approximateAsLinestring(this->numberOfControlPoints(), REQUIRED_SEGMENT_LENGTH_ACCURACY,
                                  cumulate_segment_lengths);
    return length;
  }

  /**
   *  @brief Generate an approximation of curve as a line string
   *  NOTE: This might throw std::exception or std::bad_alloc, if heap allocation fails
   *
   *  @return a vector containing vertices between line segments
   *  @return a vector containning a single vertex, if curve has only one control point
   */
  [[nodiscard]] constexpr std::vector<point> asLinestring() const {
    std::size_t constexpr INITIAL_RESERVED_SIZE = 1024u;
    real constexpr REQUIRED_SEGMENT_LENGTH_ACCURACY = real(1e-6);

    std::vector<point> segment_ends;
    segment_ends.reserve(INITIAL_RESERVED_SIZE); // Might throw
    auto const cumulate_segment_lengths = [&segment_ends](real const segment_length,
                                                          std::pair<real, point> const segment_end) -> void {
      (void)segment_length;
      segment_ends.emplace_back(segment_end.second); // Might throw
    };

    this->approximateAsLinestring(this->numberOfControlPoints(), REQUIRED_SEGMENT_LENGTH_ACCURACY,
                                  cumulate_segment_lengths);

    return segment_ends;
  }

  /**
   *  @brief Approximate the curve with linestring. Initially split curve to 'initial_vertices' vertices and then
   * continue to split segments between adjacent vertices until segment length varies at maximum of
   * 'max_segment_error'. A--B accept new vertex C,
   *                       \/
   *                       C
   *  if |AC|+|CB|>|AB|+'max_segment_error'
   *
   *    @startuml
   *    circle A
   *    circle B
   *    circle C
   *
   *    A -right- B
   *    A -- C
   *    B -- C
   *    @enduml
   *
   *  NOTE: This might throw std::exception or std::bad_alloc, if heap allocation fails
   *
   *  @param initial_vertices Amount of initial vertices before segment splitting start, minimum of 2
   *  @param max_segment_error The amount of that new added vertice must change length of segments to accept segment
   * split
   *  @param pick_segment Function or lambda function to process accepted segments Function arguments are (real)
   * segment length, segment as a pair of u value and segment end point. The segments are feed to function in
   * reverse order.
   */
  constexpr void
  approximateAsLinestring(std::size_t const initial_vertices, real const max_segment_error,
                          std::function<void(real const, std::pair<real, point> const)> pick_segment) const {
    std::vector<std::pair<real, point>> vertices;
    std::size_t constexpr INITIAL_RESERVATION = 1u << 10u;
    vertices.reserve(std::max(INITIAL_RESERVATION, initial_vertices)); // Might throw
    constexpr std::size_t MINIMUM_AMOUNT = 2u;

    // Fill initial vertices
    {
      std::size_t const N = std::max(MINIMUM_AMOUNT, initial_vertices);
      real const step = (U_MAX - U_MIN) / real(N - 1u);
      for (std::size_t i = 0u; i < N; ++i) {
        real const u = real(i) * step;
        point const p = this->C(u);
        vertices.emplace_back(u, p);
      }
    }

    // Pop segments, until no more segments available
    while (vertices.size() >= MINIMUM_AMOUNT) {
      // Cannot use reference there as allocation can invalidate reference!
      auto const segment_end = vertices.back();
      vertices.pop_back();

      // Split segment, until accurate enough
      for (;;) {
        auto const &segment_start = vertices.back();
        real const single_segment_length = segment_start.second.distance(segment_end.second);

        real const middle_u = real(0.5) * (segment_start.first + segment_end.first);

        point const middle_point = this->C(middle_u);
        real const dual_segment_length =
            segment_start.second.distance(middle_point) + middle_point.distance(segment_end.second);

        if (dual_segment_length > (single_segment_length + max_segment_error)) {
          // Try again with a new start vertex for current segment
          vertices.emplace_back(middle_u, middle_point); // Might throw
        } else {
          pick_segment(single_segment_length, segment_end);
          break;
        }
      }
    }

    // Process the last vertex
    assert(vertices.size() == 1u);
    pick_segment(U_MIN, vertices.at(0u));
  }

  /**
   *  @brief Select initial starting points for Newton-Raphson method to find a point from curve that is nearest to
   * point 'p'. @see findNearestPointFor() generates a coarse linestring for curve and calculate distances from each
   * vertex to point 'p'
   *
   *  @param[inout] nearest the candidates to use with Newton-Raphson method as a pair of squared distance and 'u'
   * Invalid guesses are marked as pair of <inf,inf>
   *  @param p The point that is used to calculate (squared) distances
   *  @param curveApproximation segment length error for generation of line string
   */
  constexpr void initialGuessesFromCurve(std::vector<std::pair<real, real>> &nearest, point const p,
                                         real const curveApproximation) const noexcept {
    for (auto &p : nearest) {
      p = std::make_pair(std::numeric_limits<real>::infinity(), std::numeric_limits<real>::infinity());
    }

    auto const pickBestInitialGuesses = [this, p, &nearest](real const segment_length,
                                                            std::pair<real, point> const segment_end) -> void {
      (void)segment_length;
      point const difference = segment_end.second - p;
      real const distance_squared = difference.lengthSquared();
      std::size_t const N = nearest.size();
      if (distance_squared < nearest.at(N - 1u).first) {
        std::size_t i = N - 1u;
        for (; i < N; --i) {
          if (distance_squared < nearest.at(i).first) {
            continue;
          }
          break;
        }
        // Now this element 'i' of array contains a value that is better
        // Keep it
        ++i;

        // Move worse results
        for (std::size_t j = N - 1u; j > i; --j) {
          nearest.at(j) = nearest.at(j - 1u);
        }

        // Save the new value
        nearest.at(i).first = distance_squared;
        nearest.at(i).second = segment_end.first;
      }
    };

    this->approximateAsLinestring(this->numberOfControlPoints(), curveApproximation, pickBestInitialGuesses);
  }

  /**
   *  @brief Find one of nearest point from curve to point 'p'.
   *  See documentation for @see closestCurvePointFor()
   *  NOTE: This might throw std::exception or std::bad_alloc, if heap allocation fails
   *
   *  @param p The point used to calculate distances to the curve
   *  @return one of the nearest point from curve to point 'p'l
   */
  [[nodiscard]] constexpr point findNearestPointFor(point const p) const {
    real constexpr REQUIRED_SEGMENT_LENGTH_ACCURACY = real(1e-3);

    // Make a grude approximation of curve as a linestring and pick the top
    // 2N-1 of vertices as starting points for Newton method
    std::size_t const N = (this->numberOfControlPoints() << 1u) - 1u;
    std::vector<std::pair<real, real>> nearest(N); // Might throw
    initialGuessesFromCurve(nearest, p, REQUIRED_SEGMENT_LENGTH_ACCURACY);

    /**
       Define distance squared function:
       \f$D_{p}(u)=\left|C(u)-p\right|^{2}=\left|C(u)-p\right|_{x}^{2}+\cdots\f$
    */
    auto const distanceSquared = [p, this](real const u) -> real {
      point const C_u = this->C(u);

      point const difference = C_u - p;
      real const distance_squared = difference.lengthSquared();
      return distance_squared;
    };

    // Define the 1st derivate:
    // f^2 = 2*f*f'
    // D_p'(u)/du = 2 * (C(u) - p) * C'(u)
    //            = 2 * (C(u).x - p.x) * C'(u).x + ...
    auto const dDistanceSquared = [p, this](real const u) -> real {
      point const C_u = this->C(u);

      point const difference = (C_u - p) * real(2);
      point const derivate = this->dC(u);
      real const result = difference * derivate;
      return result;
    };

    /**
       Define the \f$2^{nd}\f$ derivate
       Approximate from derivate defined by limit
       \f$\frac{d^{2}}{du^{2}}D_{p}(u)=\underset{h\rightarrow0}{\lim}\frac{\frac{d}{du}D_{p}(u+h)-\frac{d}{du}D_{p}(u)}{h}\approx\frac{\triangle\frac{d}{du}D_{p}(u)}{\triangle
       u}\f$
    */
    auto const d2DistanceSquared = [dDistanceSquared](real const u) -> real {
      constexpr real du = real(1e-5);
      real const u_m = std::max(U_MIN, u - du);
      real const u_p = std::min(U_MAX, u + du);
      real const Du_m = dDistanceSquared(u_m);
      real const Du_p = dDistanceSquared(u_p);
      real const result = (Du_p - Du_m) / (u_p - u_m);
      return result;
    };

    // Limit the search to valid range of u and maximum iteration rounds
    real constexpr LOWER_BOUND = U_MIN;
    real constexpr UPPER_BOUND = U_MAX;
    std::size_t constexpr NEWTON_MAX_ROUND_LIMIT = 20u;

    // Reserve space for potential locations of minimum as tuples of distance^2 and u, might throw
    std::vector<std::pair<real, real>> potential_minimums(
        N + 2u, std::make_pair(std::numeric_limits<real>::infinity(), std::numeric_limits<real>::infinity()));

    // Minimum can be at end of ranges of u = [0;1]
    potential_minimums.at(0u) = std::make_pair(distanceSquared(U_MIN), U_MIN);
    potential_minimums.at(1u) = std::make_pair(distanceSquared(U_MAX), U_MAX);

    // Minimum can be at location, where derivate of distance (squared) function is zero. Use Newton-Raphson method
    // to find them by using the 1st and 2nd derivates (approximate) of distance (squared) function
    for (std::size_t i = 0u; i < N; ++i) {
      if (std::isfinite(nearest.at(i).second)) {
        std::optional<real> const has_root = curve::internal::NewtonRaphson<real>::findRoot(
            LOWER_BOUND, nearest.at(i).second, UPPER_BOUND, NEWTON_MAX_ROUND_LIMIT, dDistanceSquared,
            d2DistanceSquared);
        if (has_root.has_value()) {
          potential_minimums.at(i + 2u) =
              std::make_pair(distanceSquared(has_root.value()), potential_minimums.at(i).second = has_root.value());
        }
      }
    }

    // Sort (asc) by distances
    std::sort(
        potential_minimums.begin(), potential_minimums.end(),
        [](std::pair<real, real> const &a, std::pair<real, real> const &b) -> bool { return a.first < b.first; });

    // There shall be at least two finite results from u=0 and u=1
    return this->C(potential_minimums.at(0u).second);
  }

private:
  span const controlPoints_;
};
} // namespace internal

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
requires std::is_same_v<CP, ControlPoint<typename CP::point>>
struct DelegationInterface {
  using ControlPoint = CP;
  using point = typename ControlPoint::point;
  using real = typename point::real;

  using point_or_issue = std::variant<point, ValidityIssue>;
  using real_or_issue = std::variant<real, ValidityIssue>;
  using vector_of_points_or_issue = std::variant<std::vector<point>, ValidityIssue>;

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
  [[nodiscard]] constexpr point_or_issue closestCurvePointFor(point const p) const noexcept;
  /// @}
};

/**
 *  @brief The class to validate rational bezier curves with control point of type 'CP', using either Point2<> or
 * Point3<> and if curve is valid, calculate points from it.
 */
template <class CP>
requires std::is_same_v<CP, ControlPoint<typename CP::point>>
class ValidateRational : private internal::CalculateRational<const CP>, public DelegationInterface<CP> {
private:
  using interface = DelegationInterface<CP>;

public:
  using ControlPoint = const CP; //< This class uses only read-only control points
  using point = typename ControlPoint::point;
  using real = typename point::real;
  using span = typename std::span<ControlPoint const>;

  using point_or_issue = typename interface::point_or_issue;
  using real_or_issue = typename interface::real_or_issue;
  using vector_of_points_or_issue = typename interface::vector_of_points_or_issue;

  /**
   *  @brief Move constructor for object to validate and calculate rational bezier curve from given control points
   *
   *  @param controlPoints read-only span from actual container of control points
   */
  inline constexpr ValidateRational(span const &&controlPoints) noexcept : calculator(std::move(controlPoints)) {}

  /**
   *  @brief A default destructor
   */
  inline constexpr ~ValidateRational() noexcept = default;

private:
  using calculator =
      internal::CalculateRational<ControlPoint>; //< Calculator class uses only read-only control points

  /**
   *  @brief Check, if curve parameter 'u' is in valid range [0;1]
   *
   *  @param u curve parameter, valid range is [0;1]
   *  @return std::nullopt, if parameter 'u' is valid
   *  @return ValidityIssue::ISSUE_U_IS_INVALID as std::optional otherwise
   */
  [[nodiscard]] static inline constexpr std::optional<ValidityIssue> isValidU(real const u) noexcept {
    return ((u >= calculator::U_MIN) && (u <= calculator::U_MAX))
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
    return (this->numberOfControlPoints() >= count)
               ? std::nullopt
               : std::make_optional(ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS);
  }

  /**
   *  @brief Check, if sum of curve's weights of control points isn't too close to zero
   *
   *  @param u curve parameter, valid range is [0;1]
   *  @return std::nullopt, if sum of weights differs from zero more than real::epsilon
   *  @return ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS as std::optional otherwise
   */
  [[nodiscard]] inline constexpr std::optional<ValidityIssue> hasValidWeights(real const u) const noexcept {
    real const sum_w = [u, this]() -> real {
      span const &cp = this->getSpan();
      std::size_t const n = cp.size() - 1u;
      real sum_w = real(0);
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::B(i, n, u) * cp[i].w();
        sum_w += w;
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
    std::optional<ValidityIssue> const issue = this->hasEnoughControlPoints(count);
    if (issue.has_value()) {
      return issue.value();
    }
    return this->hasValidWeights(real(0.5));
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
    return this->hasValidWeights(u);
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

    for (real const u : us) {
      issue = this->isValidU(u);
      if (issue.has_value()) {
        return issue.value();
      }

      issue = this->hasValidWeights(u);
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

    return this->C(u);
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

    std::vector<point> points;
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

    if (this->numberOfControlPoints() < 2u) {
      // point as a degenerate curve
      return real(0);
    } else {
      return this->length();
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

    if (this->numberOfControlPoints() < 2u) {
      // point as a degenerate curve
      std::vector<point> single_point(1u, this->getSpan().front().p());
      return single_point;
    } else {
      std::vector<point> lineSegments = this->asLinestring();
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

    return this->dC(u);
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

    point const velocity = this->dC(u);
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

    point const velocity = this->dC(u);

    // It should be impossible NOT to normalize tangent ..
    std::optional<point> const has_unit_length = velocity.normalize();
    // .. but, if that happens, return vector with zero length
    constexpr const point zero_length = point();
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
  [[nodiscard]] inline constexpr point_or_issue closestCurvePointFor(point const p) const noexcept {
    std::optional<ValidityIssue> const issue = this->isValid(1u);
    if (issue.has_value()) {
      return issue.value();
    }

    if (this->numberOfControlPoints() == 1u) {
      // point as a degenerate curve
      return this->getSpan()[0u].p();
    } else {
      return this->findNearestPointFor(p);
    }
  }
};

/**
 *  @brief The class for rational bezier curves with management of control points of type 'CP', either Point2<> or
 * Point3<>. Delegates calls to @ref curve::bezier::rational::ValidateRational
 */
template <class CP>
requires std::is_same_v<CP, ControlPoint<typename CP::point>>
class Rational : public DelegationInterface<CP> {
private:
  using interface = DelegationInterface<CP>;

public:
  using ControlPoint = typename interface::ControlPoint;
  using point = typename interface::point;
  using real = typename interface::point::real;
  using container = typename std::vector<ControlPoint>;
  using span = typename std::span<ControlPoint const>;
  using validate = ValidateRational<ControlPoint>;

  using point_or_issue = typename interface::point_or_issue;
  using real_or_issue = typename interface::real_or_issue;
  using vector_of_points_or_issue = typename interface::vector_of_points_or_issue;

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
  inline constexpr Rational(container const &cp) {
    this->controlPoints_.reserve(cp.size()); // Might throw
    this->controlPoints_.insert(this->controlPoints_.begin(), cp.cbegin(), cp.cend());
  }

  /**
   *  @brief Constructor with vector of control points moved to the curve.
   *  NOTE: cumulative weight of control points must not be near zero.
   *
   *  @param cp The control points of 'CP to be moved to this curve
   */
  inline constexpr Rational(container &&cp) noexcept { this->controlPoints_ = std::move(cp); }

  /**
   *  @brief Destructor for rational bezier curve
   */
  inline constexpr ~Rational() noexcept { this->controlPoints_ = container(); }

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
    validate const validator(this->controlPointContainerAsSpan());
    return validator.pointAt(u);
  }

  /**
   *  @brief Calculate the length of the curve
   *
   *  @return the length of curve as std::variant
   *  @return ValidityIssue as std::variant, if any problem with curve
   */
  [[nodiscard]] inline constexpr real_or_issue curveLength() const noexcept {
    validate const validator(this->controlPointContainerAsSpan());
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
    validate const validator(this->controlPointContainerAsSpan());
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
    validate const validator(this->controlPointContainerAsSpan());
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
    validate const validator(this->controlPointContainerAsSpan());
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
    validate const validator(this->controlPointContainerAsSpan());
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
  [[nodiscard]] inline point_or_issue closestCurvePointFor(point const p) const noexcept {
    validate const validator(this->controlPointContainerAsSpan());
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
  [[nodiscard]] inline constexpr span controlPointContainerAsSpan() const noexcept {
    return span(this->controlPoints_.cbegin(), this->controlPoints_.size());
  }

protected:
  container controlPoints_; //< The container for control points
};

} // namespace rational
} // namespace bezier
} // namespace curve

#endif
