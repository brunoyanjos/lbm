#pragma once

#include <cstdint>

#include "core/types.cuh"
#include "core/physics.h"

#include "lbm/domain/mask_utils.cuh"
#include "lbm/hermite/hermite.cuh"

#include "lbm/boundary/common/factor.cuh"

#include "lbm/boundary/fluid/solve_newton.cuh"
#include "lbm/boundary/fluid/accumulator.cuh"
#include "lbm/boundary/fluid/evaluate_incoming.cuh"
#include "lbm/boundary/fluid/evaluate_outgoing.cuh"
#include "lbm/boundary/fluid/eval_density.cuh"
#include "lbm/boundary/fluid/build_system.cuh"
#include "lbm/boundary/fluid/variables.cuh"

namespace boundary::fluid
{
    __device__ inline void apply_boundary(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMoments &M)
    {
        const mask_t outgoing_mask = valid_mask;
        const mask_t incoming_mask = mask_opp(valid_mask);

        Accumulator acc{};

#pragma unroll
        for (int i = 0; i < Stencil::Q; ++i)
        {
            if (dir_valid(incoming_mask, i))
                evaluate_incoming(acc, pop, i);

            if (dir_valid(outgoing_mask, i))
                evaluate_outgoing(acc, i);
        }

        acc.in.normalize();

        SystemData<UnknownMomentList::size, NonLinearSystemList::size> S{};

        build_system(S, acc);
        solve_newton(S, M);

        M.rho = eval_density(acc, M);
    }
}