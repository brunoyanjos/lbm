#pragma once

#include "properties.cuh"

#include "../../core/types.cuh"
#include "../../core/math_utils.cuh"

namespace SQUARE_CAVITY
{
    __device__ __forceinline__ void bc_velocity(int x, int y, real_t &ux, real_t &uy)
    {
        ux = r::zero;
        uy = r::zero;

        if (y == NY - 1)
        {
            ux = U_MAX;
        }
    }
}