#pragma once

#include "core/types.cuh"
#include "core/math_utils.cuh"

#include "lbm/moment/node_moments.cuh"

__device__ __forceinline__ real_t rel_step(real_t x_new, real_t x_old)
{
    const real_t eps = real_t(1e-12);
    const real_t denom = fmax(r_abs(x_new), eps);
    return r_abs(x_new - x_old) / denom;
}

__device__ __forceinline__
    real_t
    max_rel_step(const NodeMoments &a, const NodeMoments &b)
{
    real_t error = rel_step(a.ux, b.ux);
    error = fmax(error, rel_step(a.uy, b.uy));
    error = fmax(error, rel_step(a.mxx, b.mxx));
    error = fmax(error, rel_step(a.mxy, b.mxy));
    error = fmax(error, rel_step(a.myy, b.myy));
    return error;
}