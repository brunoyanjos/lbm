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
    template <int RegOrder, bool Rec, bool HighOrder>
    __device__ inline void apply_boundary(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
    {
        const mask_t outgoing_mask = valid_mask;
        const mask_t incoming_mask = mask_opp(valid_mask);

        Accumulator<RegOrder, Rec> acc{};

#pragma unroll
        for (int i = 0; i < Stencil::Q; ++i)
        {
            if (dir_valid(incoming_mask, i))
                evaluate_incoming(acc, pop, i);

            if (dir_valid(outgoing_mask, i))
                evaluate_outgoing(acc, i);
        }

        acc.in.normalize();

        SystemData<UnknownMomentList::size, FluidNonlinearMoments<RegOrder, Rec>::size> S{};

        build_system(S, acc);
        solve_newton(S, M);

        M.rho = eval_density(acc, M);
    }

    template <bool HighOrder>
    __device__ inline void apply_boundary(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMomentsFor<3, false, HighOrder> &M)
    {
        const mask_t outgoing_mask = valid_mask;
        const mask_t incoming_mask = mask_opp(valid_mask);

        Accumulator<3, false> acc{};

#pragma unroll
        for (int i = 0; i < Stencil::Q; ++i)
        {
            if (dir_valid(incoming_mask, i))
                evaluate_incoming(acc, pop, i);

            if (dir_valid(outgoing_mask, i))
                evaluate_outgoing(acc, i);
        }

        acc.in.normalize();

        SystemData<UnknownMomentList::size, FluidNonlinearMoments<3, false>::size> S{};

        build_system(S, acc);
        solve_newton(S, M);

        if constexpr (HighOrder)
        {
            M.mxxx = M.ux * M.ux * M.ux;
            M.myyy = M.uy * M.uy * M.uy;
        }

        M.mxxy = M.ux * M.ux * M.uy;
        M.mxyy = M.ux * M.uy * M.uy;

        M.rho = eval_density(acc, M);
    }
}