#pragma once

#include "../../../core/types.cuh"

struct FluidNewtonCoefficients
{
    real_t A[5 * 8]{};
    real_t b[5]{};

    __device__ __forceinline__ real_t &coeff(int r, int c)
    {
        return A[r * 8 + c];
    }

    __device__ __forceinline__ const real_t &coeff(int r, int c) const
    {
        return A[r * 8 + c];
    }
};