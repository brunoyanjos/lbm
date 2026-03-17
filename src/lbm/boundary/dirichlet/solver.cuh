#pragma once

#include <cstdint>
#include <cstdio>

#include "core/types.cuh"

#include "lbm/domain/mask_utils.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/accumulate_incoming.cuh"
#include "lbm/boundary/common/accumulate_outgoing.cuh"
#include "lbm/boundary/common/accumulator.cuh"
#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/common/gauss_elimination.cuh"
#include "lbm/boundary/common/layout/unknown_moments.cuh"
#include "lbm/boundary/common/eval_density.cuh"

#include "lbm/boundary/dirichlet/build_boundary_system.cuh"

__device__ __forceinline__ void apply_boundary(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    NodeMoments &M)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    MomentAccumulator<false> acc{};

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        if (dir_valid(incoming_mask, i))
            accumulate_incoming(acc, pop, i);

        if (dir_valid(outgoing_mask, i))
            accumulate_outgoing(acc, i);
    }

    acc.in.normalize();

    SystemData<UnknownMomentList<false>::size> S{};

    build_boundary_system(S, M.ux, M.uy, acc);

    gaussianElimination(S);

    M.mxx = S.x[0];
    M.mxy = S.x[1];
    M.myy = S.x[2];

    M.rho = eval_density(acc, M);
}