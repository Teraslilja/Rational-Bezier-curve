
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

using real = float;

using namespace curve::bezier::utilities;

struct Tests : BernsteinPolynomials<real>
{
    using Base = BernsteinPolynomials<real>;
    
  [[nodiscard]] static inline constexpr std::size_t binomialRecursiveMult(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialRecursiveMult(n, k);
  }

  [[nodiscard]] static inline std::size_t binomialRecursiveSum(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialRecursiveSum(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialFallingFactorial(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialFallingFactorial(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialNaive(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialNaive(n, k);
  }

  [[nodiscard]] static inline constexpr std::size_t binomialInPlace(std::size_t const n, std::size_t const k) noexcept {
    return Base::binomialInPlace(n, k);
  }
};


static void BM_recursiveMult( benchmark::State& state)
{
    std::size_t const n = state.range(0);
    std::size_t const k = n >> 1;

    for (auto _ : state)
    {
        for( size_t i = 0; i < SLOWDOWN; ++i )
        {
            std::size_t result = Tests::binomialRecursiveMult( n, k );
            benchmark::DoNotOptimize(result);
        }
    }
    state.SetComplexityN(n);
}

[[maybe_unused]] static void BM_recursiveSum( benchmark::State& state)
{
    std::size_t const n = state.range(0);
    std::size_t const k = n >> 1;

    for (auto _ : state)
    {
        for( size_t i = 0; i < SLOWDOWN; ++i )
        {
            std::size_t result = Tests::binomialRecursiveSum( n, k );
            benchmark::DoNotOptimize(result);
        }
    }
    state.SetComplexityN(n);
}

static void BM_falling_factorial( benchmark::State& state)
{
    std::size_t const n = state.range(0);
    std::size_t const k = n >> 1;

    for (auto _ : state)
    {
        for( size_t i = 0; i < SLOWDOWN; ++i )
        {
            std::size_t result = Tests::binomialFallingFactorial( n, k );
            benchmark::DoNotOptimize(result);
        }
    }
    state.SetComplexityN(n);
}

static void BM_naive( benchmark::State& state)
{
    std::size_t const n = state.range(0);
    std::size_t const k = n >> 1;

    for (auto _ : state)
    {
        for( size_t i = 0; i < SLOWDOWN; ++i )
        {
            std::size_t result = Tests::binomialNaive( n, k );
            benchmark::DoNotOptimize(result);
        }
    }
    state.SetComplexityN(n);
}

static void BM_in_place( benchmark::State& state)
{
    std::size_t const n = state.range(0);
    std::size_t const k = n >> 1;

    for (auto _ : state)
    {
        for( size_t i = 0; i < SLOWDOWN; ++i )
        {
            std::size_t result = Tests::binomialInPlace( n, k );
            benchmark::DoNotOptimize(result);
        }
    }
    state.SetComplexityN(n);
}


#define STATISTICS(R) Repetitions(R)->ComputeStatistics("min", [](std::vector<double> const& v) -> double { return *std::ranges::min_element(v); })->ComputeStatistics("max", [](std::vector<double> const& v) -> double { return *std::ranges::max_element(v); })

#define EXEC1(R,N) DenseRange(0, (N), 1)->STATISTICS((R))
#define EXEC2(N) DenseRange(0, (N), 1)->Complexity(benchmark::oAuto)
#define EXEC3(R,N) Range((N), (N))->Repetitions((R))

#if 0
BENCHMARK(BM_recursiveMult)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
BENCHMARK(BM_recursiveSum)->EXEC1(3, MIN_BENCH_SIZE);
BENCHMARK(BM_falling_factorial)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
BENCHMARK(BM_naive)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
BENCHMARK(BM_in_place)->EXEC1(REPEATIONS, MAX_BENCH_SIZE);
#endif

BENCHMARK(BM_recursiveMult)->EXEC2(67);
//BENCHMARK(BM_recursiveSum)->EXEC2(0);
BENCHMARK(BM_falling_factorial)->EXEC2(29);
BENCHMARK(BM_naive)->EXEC2(20);
BENCHMARK(BM_in_place)->EXEC2(67);

BENCHMARK(BM_recursiveMult)->EXEC3(9,20);
BENCHMARK(BM_falling_factorial)->EXEC3(9,20);
BENCHMARK(BM_naive)->EXEC3(9,20);
BENCHMARK(BM_in_place)->EXEC3(9,20);

BENCHMARK(BM_recursiveMult)->EXEC3(9,29);
BENCHMARK(BM_falling_factorial)->EXEC3(9,29);
BENCHMARK(BM_in_place)->EXEC3(9,29);

BENCHMARK(BM_recursiveMult)->EXEC3(9,67);
BENCHMARK(BM_in_place)->EXEC3(9,67);



BENCHMARK_MAIN();
