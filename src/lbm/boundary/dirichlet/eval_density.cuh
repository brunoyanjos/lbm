#pragma once

#include "core/types.cuh"

#include "lbm/boundary/dirichlet/accumulator.cuh"
#include "lbm/moment/node_moments.cuh"

namespace boundary::dirichlet
{
    [[nodiscard]] __device__ __forceinline__ real_t eval_density(const Accumulator &acc, const NodeMoments &M)
    {
        const real_t rho_denominator = acc.rho.constant + M.mxy * acc.rho.mxy;

        const real_t inv_rho = r::one / rho_denominator;

        return acc.in.rho * inv_rho;
    }
}