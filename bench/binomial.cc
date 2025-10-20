
#include "bernstein_polynomials.hpp"

#include <benchmark/benchmark.h>

#include <iostream>

//
// Configuration for benchmarks
//

static std::size_t constexpr MAX_BENCH_SIZE = 20;
static std::size_t constexpr MIN_BENCH_SIZE = 5;
static std::size_t constexpr REPEATIONS = 3u;
static std::size_t constexpr SLOWDOWN = 1000;

using namespace curve::bezier::utilities::binomial;

struct Tests
{
  [[nodiscard]] static inline constexpr std::size_t binomialNaive(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return naive::binomial(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialFallingFactorial(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return falling_factorial::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t
  binomialMultiplicationWithoutRecursion(std::size_t const n,
                                         std::size_t const k) noexcept
  {
    return multiplication_without_recursion::binomial(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t
  binomialMultiplicationWithRecursion(std::size_t const n,
                                      std::size_t const k) noexcept
  {
    return multiplication_with_recursion::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t binomialSumWithoutRecursion(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return sum_without_recursion::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t binomialSumWithRecursion(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return sum_with_recursion::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t binomial(std::size_t const n,
                                                   std::size_t const k) noexcept
  {
    return curve::bezier::utilities::binomial::binomial(n, k);
  }
};

static void
BM_Naive(benchmark::State& state)
{
  std::size_t const n = state.range(0);
  std::size_t const k = n >> 1;

  for (auto _ : state) {
    for (size_t i = 0; i < SLOWDOWN; ++i) {
      std::size_t result = Tests::binomialNaive(n, k);
      benchmark::DoNotOptimize(result);
    }
  }
  state.SetComplexityN(n);
}

static void
BM_FallingFactorial(benchmark::State& state)
{
  std::size_t const n = state.range(0);
  std::size_t const k = n >> 1;

  for (auto _ : state) {
    for (size_t i = 0; i < SLOWDOWN; ++i) {
      std::size_t result = Tests::binomialFallingFactorial(n, k);
      benchmark::DoNotOptimize(result);
    }
  }
  state.SetComplexityN(n);
}

static void
BM_MultiplicationWithoutRecursion(benchmark::State& state)
{
  std::size_t const n = state.range(0);
  std::size_t const k = n >> 1;

  for (auto _ : state) {
    for (size_t i = 0; i < SLOWDOWN; ++i) {
      std::size_t result = Tests::binomialMultiplicationWithoutRecursion(n, k);
      benchmark::DoNotOptimize(result);
    }
  }
  state.SetComplexityN(n);
}

static void
BM_MultiplicationWithRecursion(benchmark::State& state)
{
  std::size_t const n = state.range(0);
  std::size_t const k = n >> 1;

  for (auto _ : state) {
    for (size_t i = 0; i < SLOWDOWN; ++i) {
      std::size_t result = Tests::binomialMultiplicationWithRecursion(n, k);
      benchmark::DoNotOptimize(result);
    }
  }
  state.SetComplexityN(n);
}

static void
BM_SumWithoutRecursion(benchmark::State& state)
{
  std::size_t const n = state.range(0);
  std::size_t const k = n >> 1;

  for (auto _ : state) {
    for (size_t i = 0; i < SLOWDOWN; ++i) {
      std::size_t result = Tests::binomialSumWithoutRecursion(n, k);
      benchmark::DoNotOptimize(result);
    }
  }
  state.SetComplexityN(n);
}

static void
BM_SumWithRecursion(benchmark::State& state)
{
  std::size_t const n = state.range(0);
  std::size_t const k = n >> 1;

  for (auto _ : state) {
    for (size_t i = 0; i < SLOWDOWN; ++i) {
      std::size_t result = Tests::binomialSumWithRecursion(n, k);
      benchmark::DoNotOptimize(result);
    }
  }
  state.SetComplexityN(n);
}

static void
BM_Binomial(benchmark::State& state)
{
  std::size_t const n = state.range(0);
  std::size_t const k = n >> 1;

  for (auto _ : state) {
    for (size_t i = 0; i < SLOWDOWN; ++i) {
      std::size_t result = Tests::binomial(n, k);
      benchmark::DoNotOptimize(result);
    }
  }
  state.SetComplexityN(n);
}

#define STATISTICS(R)                                                          \
  Repetitions(R)                                                               \
    ->ComputeStatistics("min",                                                 \
                        [](std::vector<double> const& v) -> double {           \
                          return *std::ranges::min_element(v);                 \
                        })                                                     \
    ->ComputeStatistics("max", [](std::vector<double> const& v) -> double {    \
      return *std::ranges::max_element(v);                                     \
    })

#define EXEC1(R, N) DenseRange(0, (N), 1)->STATISTICS((R))
#define EXEC2(N) DenseRange(0, (N), 1)->Complexity(benchmark::oAuto)
#define EXEC3(N) DenseRange((N), (N), 1)->Repetitions(9)

#if 0
BENCHMARK(BM_Naive)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
BENCHMARK(BM_FallingFactorial)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
BENCHMARK(BM_MultWithoutRecursion)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
BENCHMARK(BM_MultWithRecursion)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
BENCHMARK(BM_SumWithoutRecursion)->EXEC1(3, MIN_BENCH_SIZE);
BENCHMARK(BM_SumWithRecursion)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
BENCHMARK(BM_Binomial)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
#endif

BENCHMARK(BM_Naive)->EXEC2(naive::maximumAllowedN);
BENCHMARK(BM_FallingFactorial)->EXEC2(falling_factorial::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithoutRecursion)
  ->EXEC2(multiplication_without_recursion::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithRecursion)
  ->EXEC2(multiplication_with_recursion::maximumAllowedN);
BENCHMARK(BM_SumWithoutRecursion)
  ->EXEC2(sum_without_recursion::maximumAllowedN);
BENCHMARK(BM_SumWithRecursion)
  ->EXEC2(std::min(std::size_t(25), sum_with_recursion::maximumAllowedN));
BENCHMARK(BM_Binomial)->EXEC2(sum_without_recursion::maximumAllowedN);

BENCHMARK(BM_Naive)->EXEC3(5);
BENCHMARK(BM_FallingFactorial)->EXEC3(5);
BENCHMARK(BM_MultiplicationWithoutRecursion)->EXEC3(5);
BENCHMARK(BM_MultiplicationWithRecursion)->EXEC3(5);
BENCHMARK(BM_SumWithoutRecursion)->EXEC3(5);
BENCHMARK(BM_SumWithRecursion)->EXEC3(5);
BENCHMARK(BM_Binomial)->EXEC3(5);

BENCHMARK(BM_Naive)->EXEC3(20);
BENCHMARK(BM_FallingFactorial)->EXEC3(20);
BENCHMARK(BM_MultiplicationWithoutRecursion)->EXEC3(20);
BENCHMARK(BM_MultiplicationWithRecursion)->EXEC3(20);
BENCHMARK(BM_SumWithoutRecursion)->EXEC3(20);
// BENCHMARK(BM_SumWithRecursion)->EXEC3(20); Too slow
BENCHMARK(BM_Binomial)->EXEC3(20);

// BENCHMARK(BM_Naive)->EXEC3(29);    Too large N
BENCHMARK(BM_FallingFactorial)->EXEC3(29);
BENCHMARK(BM_MultiplicationWithoutRecursion)->EXEC3(29);
BENCHMARK(BM_MultiplicationWithRecursion)->EXEC3(29);
BENCHMARK(BM_SumWithoutRecursion)->EXEC3(29);
// BENCHMARK(BM_SumWithRecursion)->EXEC3(29);
BENCHMARK(BM_Binomial)->EXEC3(29);

// BENCHMARK(BM_Naive)->EXEC3(67);
// BENCHMARK(BM_FallingFactorial)->EXEC3( 67);
BENCHMARK(BM_MultiplicationWithoutRecursion)->EXEC3(67);
BENCHMARK(BM_MultiplicationWithRecursion)->EXEC3(67);
BENCHMARK(BM_SumWithoutRecursion)->EXEC3(67);
// BENCHMARK(BM_SumWithRecursion)->EXEC3(67);
BENCHMARK(BM_Binomial)->EXEC3(67);

BENCHMARK_MAIN();
