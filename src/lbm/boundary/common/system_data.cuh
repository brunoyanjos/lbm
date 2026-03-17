#pragma once

#include <cstddef>

template <std::size_t N, std::size_t E = 0>
struct SystemData
{
    static constexpr std::size_t rows = N;
    static constexpr std::size_t cols = N + E;

    real_t A[rows * cols]{};
    real_t b[rows]{};
    real_t x[rows]{};

    __device__ __forceinline__ real_t &coeff(int i, int j)
    {
        return A[i * cols + j];
    }

    __device__ __forceinline__ const real_t &coeff(int i, int j) const
    {
        return A[i * cols + j];
    }
};