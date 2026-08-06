#pragma once
#include <type_traits>
#include <cmath>
#include "types.cuh"

template <typename T>
__host__ __device__ __forceinline__ T dsqrt(T x)
{
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "dsqrt only supports float or double");
#if defined(__CUDA_ARCH__)
    if constexpr (std::is_same_v<T, float>)
        return ::sqrtf(x);
    else
        return ::sqrt(x);
#else
    using std::sqrt;
    return sqrt(x);
#endif
}

template <typename T>
__host__ __device__ __forceinline__ constexpr T constexpr_sqrt(T x)
{
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "constexpr_sqrt only supports float or double");

    if (x <= T{0})
        return T{0};

    T guess = x;
    T previous = T{0};
    for (int i = 0; i < 64 && guess != previous; ++i)
    {
        previous = guess;
        guess = (guess + x / guess) * T{0.5};
    }
    return guess;
}

template <typename T>
__host__ __device__ __forceinline__ constexpr T dabs(T x)
{
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "dabs only supports float or double");
    return x < T{0} ? -x : x;
}

// Conveniência: versões específicas pra real_t
__host__ __device__ __forceinline__ real_t r_sqrt(real_t x) { return dsqrt<real_t>(x); }
__host__ __device__ __forceinline__ real_t r_abs(real_t x) { return dabs<real_t>(x); }
