//
// (C) Matti Lehtonen 2023
//

#include "newton_raphson.hpp"

#include <gtest/gtest.h>

#include <iostream>

using real = float;
using NewtonRaphson = curve::internal::NewtonRaphson<real>;

struct NewtonRaphsonData {
  std::function<std::optional<real>(real const)> const f;
  std::function<std::optional<real>(real const)> const df;
  real const lowerBound;
  real const upperBound;
  std::size_t const roundLimit;
  real const initialU;

  std::optional<real> const expectedResult;

  inline friend std::ostream &operator<<(std::ostream &out, NewtonRaphsonData const &data) {
    out << "{";
    out << data.lowerBound << "," << data.upperBound << "," << data.roundLimit << "," << data.initialU << ",";
    if (data.expectedResult.has_value()) {
      out << data.expectedResult.value();
    } else {
      out << "'no value'";
    }
    out << "}";
    return out;
  }
};

class NewtonRaphsonTests : public ::testing::TestWithParam<NewtonRaphsonData> {};

TEST_P(NewtonRaphsonTests, findroot_correctResult) {
  static constexpr real EPS = real(1e-5);
  auto const &params = GetParam();

  auto const result = NewtonRaphson::findRoot(params.lowerBound, params.initialU, params.upperBound,
                                              params.roundLimit, params.f, params.df);

  ASSERT_EQ(result.has_value(), params.expectedResult.has_value());
  if (params.expectedResult.has_value()) {
    EXPECT_NEAR(result.value(), params.expectedResult.value(), EPS);
  }
}

auto const parabel = [](real const u) -> std::optional<real> { return u * u - real(1); };
auto const d_parabel = [](real const u) -> std::optional<real> { return real(2) * u; };

auto const constant0 = [](real const u) -> std::optional<real> {
  (void)u;
  return real(0);
};
auto const constant1 = [](real const u) -> std::optional<real> {
  (void)u;
  return real(1);
};
auto const d_constant = [](real const u) -> std::optional<real> {
  (void)u;
  return real(0);
};

auto const failure_func = [](real const u) -> std::optional<real> {
  (void)u;
  return std::nullopt;
};
auto const d_failure_func = [](real const u) -> std::optional<real> {
  (void)u;
  return std::nullopt;
};

static const NewtonRaphsonData newtonRaphsonData[]{
    {parabel, d_parabel, real(-100), real(100), 1u, real(10), std::nullopt},

    {parabel, d_parabel, real(-100), real(100), 10u, real(10), real(1)},
    {parabel, d_parabel, real(-100), real(100), 10u, real(-10), real(-1)},

    {constant1, d_constant, real(-100), real(100), 10u, real(0), real(0)},
    {constant0, d_constant, real(-100), real(100), 10u, real(0), real(0)},

    {parabel, d_failure_func, real(-100), real(100), 10u, real(10), std::nullopt},
    {failure_func, d_failure_func, real(-100), real(100), 10u, real(-10), std::nullopt},
    {failure_func, d_parabel, real(-100), real(100), 10u, real(-10), std::nullopt},

    {parabel, d_parabel, real(2), real(10), 10u, real(9), std::nullopt},
};

INSTANTIATE_TEST_SUITE_P(Fixture, NewtonRaphsonTests, ::testing::ValuesIn(newtonRaphsonData));
