#pragma once
#include "../../domain_tags.cuh"
#include "../../../../core/types.cuh"
#include "../../../../core/geometries/square_cavity/geometry.h"
#include "../../../../core/geometries/square_cavity/physics.h"
#include "../../../../core/math_utils.cuh"

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