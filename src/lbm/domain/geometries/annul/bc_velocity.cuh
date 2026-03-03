#pragma once
#include "../../domain_tags.cuh"
#include "../../../../core/types.cuh"
#include "../../../../core/physics.h"
#include "../../../../core/geometry.h"
#include "../../../../core/math_utils.cuh"

namespace ANNUL
{
    __device__ __forceinline__ void bc_velocity(int x, int y, real_t &ux, real_t &uy)
    {
        ux = r::zero;
        uy = r::zero;

        real_t c, s, r;
        polar_unit_vectors(x, y, c, s, r);

        const real_t tol = real_t(1.5);
        const real_t dr = r_abs(r - R_IN);

        if (dr <= tol)
        {
            const real_t utheta = U_WALL;
            ux = -utheta * s;
            uy = utheta * c;
        }
    }
}