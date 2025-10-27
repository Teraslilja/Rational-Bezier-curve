
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
    return Naive::binomial(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialFallingFactorial(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return FallingFactorial::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t
  binomialMultiplicationWithoutRecursion(std::size_t const n,
                                         std::size_t const k) noexcept
  {
    return MultiplicationWithoutRecursion::binomial(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t
  binomialMultiplicationWithRecursion(std::size_t const n,
                                      std::size_t const k) noexcept
  {
    return MultiplicationWithRecursion::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t binomialSumWithoutRecursion(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return SumWithoutRecursion::binomial(n, k);
  }

  [[nodiscard]] static inline std::size_t binomialSumWithRecursion(
    std::size_t const n,
    std::size_t const k) noexcept
  {
    return SumWithRecursion::binomial(n, k);
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

BENCHMARK(BM_Naive)->EXEC2(Naive::maximumAllowedN);
BENCHMARK(BM_FallingFactorial)->EXEC2(FallingFactorial::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithoutRecursion)
  ->EXEC2(MultiplicationWithoutRecursion::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithRecursion)
  ->EXEC2(MultiplicationWithRecursion::maximumAllowedN);
BENCHMARK(BM_SumWithoutRecursion)->EXEC2(SumWithoutRecursion::maximumAllowedN);
BENCHMARK(BM_SumWithRecursion)->EXEC2(SumWithRecursion::maximumAllowedN);
BENCHMARK(BM_Binomial)->EXEC2(SumWithoutRecursion::maximumAllowedN);

BENCHMARK(BM_Naive)->EXEC3(5);
BENCHMARK(BM_FallingFactorial)->EXEC3(5);
BENCHMARK(BM_MultiplicationWithoutRecursion)->EXEC3(5);
BENCHMARK(BM_MultiplicationWithRecursion)->EXEC3(5);
BENCHMARK(BM_SumWithoutRecursion)->EXEC3(5);
BENCHMARK(BM_SumWithRecursion)->EXEC3(5);
BENCHMARK(BM_Binomial)->EXEC3(5);

BENCHMARK(BM_Naive)->EXEC3(Naive::maximumAllowedN);
BENCHMARK(BM_FallingFactorial)->EXEC3(Naive::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithoutRecursion)->EXEC3(Naive::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithRecursion)->EXEC3(Naive::maximumAllowedN);
BENCHMARK(BM_SumWithoutRecursion)->EXEC3(Naive::maximumAllowedN);
BENCHMARK(BM_SumWithRecursion)->EXEC3(Naive::maximumAllowedN);
BENCHMARK(BM_Binomial)->EXEC3(Naive::maximumAllowedN);

// BENCHMARK(BM_Naive)->EXEC3(FallingFactorial::maximumAllowedN);    Too large N
BENCHMARK(BM_FallingFactorial)->EXEC3(FallingFactorial::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithoutRecursion)
  ->EXEC3(FallingFactorial::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithRecursion)
  ->EXEC3(FallingFactorial::maximumAllowedN);
BENCHMARK(BM_SumWithoutRecursion)->EXEC3(FallingFactorial::maximumAllowedN);
// BENCHMARK(BM_SumWithRecursion)->EXEC3(FallingFactorial::maximumAllowedN);
BENCHMARK(BM_Binomial)->EXEC3(FallingFactorial::maximumAllowedN);

// BENCHMARK(BM_Naive)->EXEC3(MultiplicationWithoutRecursion::maximumAllowedN);
// BENCHMARK(BM_FallingFactorial)->EXEC3(MultiplicationWithoutRecursion::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithoutRecursion)
  ->EXEC3(MultiplicationWithoutRecursion::maximumAllowedN);
BENCHMARK(BM_MultiplicationWithRecursion)
  ->EXEC3(MultiplicationWithoutRecursion::maximumAllowedN);
BENCHMARK(BM_SumWithoutRecursion)
  ->EXEC3(MultiplicationWithoutRecursion::maximumAllowedN);
// BENCHMARK(BM_SumWithRecursion)->EXEC3(MultiplicationWithoutRecursion::maximumAllowedN);
BENCHMARK(BM_Binomial)->EXEC3(MultiplicationWithoutRecursion::maximumAllowedN);

// BENCHMARK(BM_Naive)->EXEC3(SumWithoutRecursion::maximumAllowedN);
// BENCHMARK(BM_FallingFactorial)->EXEC3(SumWithoutRecursion::maximumAllowedN);
// BENCHMARK(BM_MultiplicationWithoutRecursion)->EXEC3(SumWithoutRecursion::maximumAllowedN);
// BENCHMARK(BM_MultiplicationWithRecursion)->EXEC3(SumWithoutRecursion::maximumAllowedN);
BENCHMARK(BM_SumWithoutRecursion)->EXEC3(SumWithoutRecursion::maximumAllowedN);
// BENCHMARK(BM_SumWithRecursion)->EXEC3(SumWithoutRecursion::maximumAllowedN);
BENCHMARK(BM_Binomial)->EXEC3(SumWithoutRecursion::maximumAllowedN);

BENCHMARK_MAIN();
