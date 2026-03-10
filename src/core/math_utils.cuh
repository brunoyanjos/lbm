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
__host__ __device__ __forceinline__ T dabs(T x)
{
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                  "dabs only supports float or double");
#if defined(__CUDA_ARCH__)
    if constexpr (std::is_same_v<T, float>)
        return ::fabsf(x);
    else
        return ::fabs(x);
#else
    using std::abs;
    return abs(x);
#endif
}

__host__ __device__ __forceinline__ real_t r_sqrt(real_t x) { return dsqrt<real_t>(x); }
__host__ __device__ __forceinline__ real_t r_abs(real_t x) { return dabs<real_t>(x); }

__host__ __device__ __forceinline__ void polar_unit_vectors(int x, int y,
                                                            real_t xc, real_t yc,
                                                            real_t &c, real_t &s,
                                                            real_t &r)
{
    const real_t dx = r_cast(x) - xc;
    const real_t dy = r_cast(y) - yc;

    const real_t r2 = dx * dx + dy * dy;
#if defined(__CUDA_ARCH__)
    const real_t invr = rsqrt(r2 + real_t(1e-30));
#else
    const real_t invr = real_t(1) / r_sqrt(r2 + real_t(1e-30));
#endif

    c = dx * invr;
    s = dy * invr;
    r = r2 * invr;
}