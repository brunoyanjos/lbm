#pragma once
#include "../../domain_tags.cuh"
#include "../../../../core/types.cuh"
#include "../../../../core/geometries/couette/geometry.h"
#include "../../../../core/geometries/couette/physics.h"
#include "../../../../core/math_utils.cuh"

namespace COUETTE
{
    __device__ __forceinline__ void bc_velocity(int x, int y, real_t &ux, real_t &uy)
    {
        (void)x;

        ux = r::zero;
        uy = r::zero;

        // Parede inferior fixa
        if (y == 0)
        {
            ux = r::zero;
            uy = r::zero;
            return;
        }

        // Parede superior móvel
        if (y == NY - 1)
        {
            ux = U_MAX;
            uy = r::zero;
            return;
        }
    }
}