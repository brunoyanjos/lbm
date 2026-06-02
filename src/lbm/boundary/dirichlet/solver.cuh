#pragma once

#include <cstdint>
#include <cstdio>

#include "core/types.cuh"

#include "lbm/domain/mask_utils.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/common/gauss_elimination.cuh"

#include "lbm/boundary/dirichlet/accumulator.cuh"
#include "lbm/boundary/dirichlet/eval_density.cuh"
#include "lbm/boundary/dirichlet/evaluate_incoming.cuh"
#include "lbm/boundary/dirichlet/evaluate_outgoing.cuh"

namespace boundary::dirichlet
{

    template <int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ void apply_boundary(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
    {
        const mask_t outgoing_mask = valid_mask;
        const mask_t incoming_mask = mask_opp(valid_mask);

        Accumulator acc{};

#pragma unroll
        for (int i = 0; i < Stencil::Q; ++i)
        {
            if (dir_valid(incoming_mask, i))
                evaluate_incoming(acc, M, pop, i);

            if (dir_valid(outgoing_mask, i))
                evaluate_outgoing(acc, M, i);
        }

        acc.in.normalize();

        M.mxy = (acc.rho.constant * acc.in.mxy - acc.mxy.constant) / (acc.mxy.mxy - acc.rho.mxy * acc.in.mxy);

        M.mxx = M.ux * M.ux;
        M.myy = M.uy * M.uy;

        M.rho = eval_density(acc, M);
    }

    template <bool HighOrder>
    __device__ __forceinline__ void apply_boundary(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMomentsFor<3, false, HighOrder> &M)
    {
        const mask_t outgoing_mask = valid_mask;
        const mask_t incoming_mask = mask_opp(valid_mask);

        Accumulator acc{};

#pragma unroll
        for (int i = 0; i < Stencil::Q; ++i)
        {
            if (dir_valid(incoming_mask, i))
                evaluate_incoming(acc, M, pop, i);

            if (dir_valid(outgoing_mask, i))
                evaluate_outgoing(acc, M, i);
        }

        acc.in.normalize();

        M.mxy = (acc.rho.constant * acc.in.mxy - acc.mxy.constant) / (acc.mxy.mxy - acc.rho.mxy * acc.in.mxy);

        M.mxx = M.ux * M.ux;
        M.myy = M.uy * M.uy;
        M.mxxy = M.ux * M.ux * M.uy;
        M.mxyy = M.ux * M.uy * M.uy;

        if constexpr (HighOrder)
        {
            M.mxxx = M.ux * M.ux * M.ux;
            M.myyy = M.uy * M.uy * M.uy;
        }

        M.rho = eval_density(acc, M);
    }
}