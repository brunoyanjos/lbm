#pragma once

#include "core/types.cuh"

#include "lbm/boundary/fluid/accumulator.cuh"
#include "lbm/moment/node_moments.cuh"

namespace boundary::fluid
{
    [[nodiscard]] __device__ __forceinline__ real_t eval_density(const Accumulator &acc, const NodeMoments &M)
    {
        const real_t rho_denominator = acc.rho.rho +
                                       M.ux * acc.rho.ux + M.uy * acc.rho.uy +
                                       M.ux * M.ux * acc.rho.uxux + M.ux * M.uy * acc.rho.uxuy + M.uy * M.uy * acc.rho.uyuy +
                                       M.mxx * acc.rho.mxx + M.mxy * acc.rho.mxy + M.myy * acc.rho.myy;

        const real_t inv_rho = r::one / rho_denominator;

        return acc.in.rho * inv_rho;
    }
}