#pragma once
#include "properties.cuh"

#include "../../lbm/domain/domain_tags.cuh"
#include "../../core/types.cuh"

#include "../../core/math_utils.cuh"

namespace POISEUILLE
{
    __device__ __forceinline__ void bc_velocity(int x, int y, real_t &ux, real_t &uy)
    {
        ux = r::zero;
        uy = r::zero;

        if (y == 0 || y == NY - 1)
        {
            ux = r::zero;
            uy = r::zero;
            return;
        }

        if (x == 0)
        {
            ux = U_MAX;
            uy = r::zero;
            return;
        }

        if (x == NX - 1)
        {
            ux = r::zero;
            uy = r::zero;
            return;
        }
    }
}