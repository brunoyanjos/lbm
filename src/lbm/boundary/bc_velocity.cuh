#pragma once

#include "core/types.cuh"
#include "core/geometry.h"
#include "core/physics.h"

#include "lbm/moment/node_moments.cuh"

__host__ __device__ __forceinline__ void polar_unit_vectors(int x, int y, real_t &c, real_t &s,
                                                            real_t &r)
{
    const real_t dx = r_cast(x) - xc;
    const real_t dy = r_cast(y) - yc;

    const real_t r2 = dx * dx + dy * dy;
    const real_t r_safe = r_sqrt(r2 + real_t(1e-30));
    const real_t invr = r::one / r_safe;

    c = dx * invr;
    s = dy * invr;
    r = r_safe;
}

__device__ __forceinline__ void bc_velocity(int x, int y, NodeMoments &M)
{
    M.ux = r::zero;
    M.uy = r::zero;

    // real_t c, s, r;
    // polar_unit_vectors(x, y, c, s, r);

    // const real_t tol = (R_OUT + R_IN) / 2;

    // if (r < tol)
    // {
    //     const real_t utheta = U_MAX;
    //     M.ux = -utheta * s;
    //     M.uy = utheta * c;
    // }
}