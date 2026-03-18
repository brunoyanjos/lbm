#pragma once

#include <cstdint>

#include "core/types.cuh"
#include "core/physics.h"

#include "lbm/domain/mask_utils.cuh"
#include "lbm/hermite/hermite.cuh"

#include "lbm/boundary/common/accumulate_incoming.cuh"
#include "lbm/boundary/common/accumulate_outgoing.cuh"
#include "lbm/boundary/common/accumulator.cuh"
#include "lbm/boundary/common/factor.cuh"
#include "lbm/boundary/common/layout/unknown_moments.cuh"
#include "lbm/boundary/common/layout/nonlinear_moments.cuh"
#include "lbm/boundary/common/eval_density.cuh"

#include "lbm/boundary/fluid/build_fluid_system.cuh"
#include "lbm/boundary/fluid/solve_fluid_newton.cuh"

__device__ inline void evaluate_fluid_node(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    NodeMoments &M)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    MomentAccumulator<true> acc{};

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        if (dir_valid(incoming_mask, i))
            accumulate_incoming(acc, pop, i);

        if (dir_valid(outgoing_mask, i))
            accumulate_outgoing(acc, i);
    }

    acc.in.normalize();

    SystemData<UnknownMomentList<true>::size, NonLinearSystemList::size> S{};

    build_fluid_system(S, acc);
    solve_fluid_newton(S, M);

    M.rho = eval_density(acc, M);
}
