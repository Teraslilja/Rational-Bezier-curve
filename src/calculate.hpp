//
// (C) Matti Lehtonen 2023
//

/**
 *  @file rational_bezier.hpp This file contains implementation of rational bezier class.
 */

#ifndef CALCULATE_H
#define CALCULATE_H

#include "bernstein_polynomials.hpp"
#include "calculation_interface.hpp"
#include "newton_raphson.hpp"

namespace curve {
namespace bezier {
namespace rational {
namespace {

/**
 *  @brief The class to calculate rational bezier curves with control point of type 'CP', using either Point2<> or
 * Point3<>
 */
template <class CP>
requires std::is_same_v<CP, ControlPoint<typename CP::Point>>
class CalculateRational : protected CalculationInterface<CP> {
public:
  using Interface = CalculationInterface<CP>;

  using ControlPoint = typename Interface::ControlPoint;
  using ConstControlPoint = ControlPoint const;
  using Point = typename Interface::ControlPoint::Point;
  using ConstPoint = Point const;
  using real = typename Interface::ControlPoint::Point::real;

  using ControlPointSpan = typename std::span<ConstControlPoint>;

  using PointAtCurve = typename Interface::PointAtCurve;
  using DistanceFromLocation = typename Interface::DistanceFromLocation;

  static real constexpr U_MIN = real(0); //< Minimum valid value of u
  static real constexpr U_MAX = real(1); //< Maximum valid value of u

  /**
   *  @brief Move constructor for object to calculate rational bezier curve from given control points
   *
   *  @param controlPoints read-only span from actual container of control points
   */
  inline constexpr CalculateRational(ControlPointSpan const &&controlPoints) noexcept
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
  [[nodiscard]] constexpr Point C(real const u) const noexcept {
    /**
       \f$C(u)=\frac{f(u)}{g(u)}\f$, where
       \f$f(u)=\underset{i=0}{\overset{n}{\sum}}B_{i,n}(u)w_{i}\mathbf{P}_{i}\f$
       and \f$g(u)=\underset{i=0}{\overset{n}{\sum}}B_{i,n}(u)w_{i}\f$
    */

    // Calculate f(u) and g(u)
    auto const [f_u, g_u] = [u](ControlPointSpan const cp) noexcept -> std::tuple<ConstPoint, real> {
      std::size_t const n = cp.size() - 1u;
      Point sum_wp;
      real sum_w = real(0);
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::B(i, n, u) * cp[i].w();
        sum_wp += (cp[i].p() * w);
        sum_w += w;
      }
      return {sum_wp, sum_w};
    }(this->controlPoints_);

    std::optional<ConstPoint> const C_u = f_u / g_u;
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
  [[nodiscard]] inline constexpr ControlPointSpan const &getSpan() const noexcept { return this->controlPoints_; }

  /**
   *  @brief Function C'(u) to calculate a derivate at 'u'. @see velocityAt()
   *
   *  @param u a position from curve, range [0;1]
   *  @return a point from curve, depending type of 'CP', either 2D or 3D
   */
  [[nodiscard]] constexpr Point dC(real const u) const noexcept {
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
    auto const [f_u, g_u] = [u](ControlPointSpan const &cp) noexcept -> std::tuple<ConstPoint, real> {
      std::size_t const n = cp.size() - 1u;
      real sum_w = real(0);
      Point sum_wp;
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::B(i, n, u) * cp[i].w();
        sum_w += w;
        sum_wp += (cp[i].p() * w);
      }
      return {sum_wp, sum_w};
    }(this->controlPoints_);

    // Calculate f'(u) and g'(u)
    auto const [df_u, dg_u] = [u](ControlPointSpan const &cp) noexcept -> std::tuple<ConstPoint, real> {
      std::size_t const n = cp.size() - 1u;
      Point sum_wp;
      real sum_w = real(0);
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::dB(i, n, u) * cp[i].w();
        sum_wp += cp[i].p() * w;
        sum_w += w;
      }
      return {sum_wp, sum_w};
    }(this->controlPoints_);

    // C'(u) = (f'(u) * g(u) - f(u) * g'(u)) / g(u)^2
    std::optional<ConstPoint> const dC_u = ((df_u * g_u) - (f_u * dg_u)) / (g_u * g_u);
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
  [[nodiscard]] constexpr Point d2C(ControlPointSpan const controlPoints, real const u) noexcept {
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
    auto const [f_u, g_u] = [u](ControlPointSpan const &cp) noexcept -> std::tuple<ConstPoint, real> {
      std::size_t const n = cp.size() - 1u;
      real sum_w = real(0);
      Point sum_wp;
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::B(i, n, u) * cp[i].w();
        sum_w += w;
        sum_wp += (cp[i].p() * w);
      }
      return {sum_wp, sum_w};
    }(controlPoints);

    // Calculate f'(u) and g'(u)
    auto const [df_u, dg_u] = [u](ControlPointSpan const &cp) noexcept -> std::tuple<ConstPoint, real> {
      std::size_t const n = cp.size() - 1u;
      Point sum_wp;
      real sum_w = real(0);
      for (std::size_t i = 0u; i <= n; ++i) {
        real const w = curve::bezier::internal::BernsteinPolynomials<real>::dB(i, n, u) * cp[i].w();
        sum_wp += cp[i].p() * w;
        sum_w += w;
      }
      return {sum_wp, sum_w};
    }(controlPoints);

    // Calculate f"(u) and g"(u)
    auto const [d2f_u, d2g_u] = [u](ControlPointSpan const &cp) noexcept -> std::tuple<ConstPoint, real> {
      std::size_t const n = cp.size() - 1u;
      Point sum_wp;
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
    std::optional<ConstPoint> const d2C_u =
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
                                                    PointAtCurve const segment_end) noexcept -> void {
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
  [[nodiscard]] constexpr std::vector<Point> asLinestring() const {
    std::size_t constexpr INITIAL_RESERVED_SIZE = 1024u;
    real constexpr REQUIRED_SEGMENT_LENGTH_ACCURACY = real(1e-6);

    std::vector<Point> segment_ends;
    segment_ends.reserve(INITIAL_RESERVED_SIZE); // Might throw
    auto const cumulate_segment_lengths = [&segment_ends](real const segment_length,
                                                          PointAtCurve const segment_end) noexcept -> void {
      (void)segment_length;
      segment_ends.emplace_back(segment_end.point); // Might throw
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
  constexpr void approximateAsLinestring(std::size_t const initial_vertices, real const max_segment_error,
                                         std::function<void(real const, PointAtCurve const)> pick_segment) const {
    std::vector<PointAtCurve> vertices;
    std::size_t constexpr INITIAL_RESERVATION = 1u << 10u;
    vertices.reserve(std::max(INITIAL_RESERVATION, initial_vertices)); // Might throw
    constexpr std::size_t MINIMUM_AMOUNT = 2u;

    // Fill initial vertices
    {
      std::size_t const N = std::max(MINIMUM_AMOUNT, initial_vertices);
      real const step = (U_MAX - U_MIN) / real(N - 1u);
      for (std::size_t i = 0u; i < N; ++i) {
        real const u = real(i) * step;
        ConstPoint p = this->C(u);
        vertices.emplace_back(u, p);
      }
    }

    // Pop segments, until no more segments available
    while (vertices.size() >= MINIMUM_AMOUNT) {
      // Cannot use reference there as allocation can invalidate reference!
      PointAtCurve const segment_end = vertices.back();
      vertices.pop_back();

      // Split segment, until accurate enough
      for (;;) {
        PointAtCurve const &segment_start = vertices.back();
        real const single_segment_length = segment_start.point.distance(segment_end.point);

        real const middle_u = real(0.5) * (segment_start.u + segment_end.u);

        ConstPoint middle_point = this->C(middle_u);
        real const dual_segment_length =
            segment_start.point.distance(middle_point) + middle_point.distance(segment_end.point);

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
  constexpr void initialGuessesFromCurve(std::vector<DistanceFromLocation> &nearest, ConstPoint p,
                                         real const curveApproximation) const noexcept {
    for (auto &p : nearest) {
      p = DistanceFromLocation();
    }

    auto const pickBestInitialGuesses = [p, &nearest](real const segment_length,
                                                      PointAtCurve const segment_end) noexcept -> void {
      (void)segment_length;
      ConstPoint difference = segment_end.point - p;
      real const distance_squared = difference.lengthSquared();
      std::size_t const N = nearest.size();
      if (distance_squared < nearest.at(N - 1u).distanceSquared) {
        std::size_t i = N - 1u;
        for (; i < N; --i) {
          if (distance_squared < nearest.at(i).distanceSquared) {
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
        nearest.at(i) = {distance_squared, segment_end.u};
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
  [[nodiscard]] constexpr Point findNearestPointFor(ConstPoint p) const {
    real constexpr REQUIRED_SEGMENT_LENGTH_ACCURACY = real(1e-3);

    // Make a grude approximation of curve as a linestring and pick the top
    // 2N-1 of vertices as starting points for Newton method
    std::size_t const N = (this->numberOfControlPoints() << 1u) - 1u;
    std::vector<DistanceFromLocation> nearest(N); // Might throw
    initialGuessesFromCurve(nearest, p, REQUIRED_SEGMENT_LENGTH_ACCURACY);

    /**
       Define distance squared function:
       \f$D_{p}(u)=\left|C(u)-p\right|^{2}=\left|C(u)-p\right|_{x}^{2}+\cdots\f$
    */
    auto const distanceSquared = [p, this](real const u) noexcept -> real {
      ConstPoint C_u = this->C(u);

      ConstPoint difference = C_u - p;
      real const distance_squared = difference.lengthSquared();
      return distance_squared;
    };

    // Define the 1st derivate:
    // f^2 = 2*f*f'
    // D_p'(u)/du = 2 * (C(u) - p) * C'(u)
    //            = 2 * (C(u).x - p.x) * C'(u).x + ...
    auto const dDistanceSquared = [p, this](real const u) noexcept -> real {
      ConstPoint C_u = this->C(u);

      ConstPoint difference = (C_u - p) * real(2);
      ConstPoint derivate = this->dC(u);
      real const result = difference * derivate;
      return result;
    };

    /**
       Define the \f$2^{nd}\f$ derivate
       Approximate from derivate defined by limit
       \f$\frac{d^{2}}{du^{2}}D_{p}(u)=\underset{h\rightarrow0}{\lim}\frac{\frac{d}{du}D_{p}(u+h)-\frac{d}{du}D_{p}(u)}{h}\approx\frac{\triangle\frac{d}{du}D_{p}(u)}{\triangle
       u}\f$
    */
    auto const d2DistanceSquared = [dDistanceSquared](real const u) noexcept -> real {
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
    std::vector<DistanceFromLocation> potential_minimums(N + 2u);

    // Minimum can be at end of ranges of u = [0;1]
    potential_minimums.at(0u) = {distanceSquared(U_MIN), U_MIN};
    potential_minimums.at(1u) = {distanceSquared(U_MAX), U_MAX};

    // Minimum can be at location, where derivate of distance (squared) function is zero. Use Newton-Raphson method
    // to find them by using the 1st and 2nd derivates (approximate) of distance (squared) function
    for (std::size_t i = 0u; i < N; ++i) {
      if (nearest.at(i).hasValidLocation()) {
        std::optional<real> const has_root = curve::internal::NewtonRaphson<real>::findRoot(
            LOWER_BOUND, nearest.at(i).u, UPPER_BOUND, NEWTON_MAX_ROUND_LIMIT, dDistanceSquared, d2DistanceSquared);
        if (has_root.has_value()) {
          potential_minimums.at(i + 2u) = {distanceSquared(has_root.value()),
                                           potential_minimums.at(i).u = has_root.value()};
        }
      }
    }

    // Sort (asc) by distances
    std::sort(potential_minimums.begin(), potential_minimums.end(),
              [](DistanceFromLocation const &a, DistanceFromLocation const &b) noexcept -> bool {
                return a.distanceSquared < b.distanceSquared;
              });

    // There shall be at least two finite results from u=0 and u=1
    return this->C(potential_minimums.at(0u).u);
  }

private:
  ControlPointSpan const controlPoints_;
};

} // namespace
} // namespace rational
} // namespace bezier
} // namespace curve

#endif
