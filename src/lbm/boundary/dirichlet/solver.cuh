#pragma once

#include <cstdint>
#include <cstdio>

#include "core/types.cuh"

#include "lbm/domain/mask_utils.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/accumulator.cuh"
#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/common/gauss_elimination.cuh"
#include "lbm/boundary/common/layout/unknown_moments.cuh"
#include "lbm/boundary/common/eval_density.cuh"

#include "lbm/boundary/dirichlet/build_boundary_system.cuh"
#include "lbm/boundary/dirichlet/dirichlet_accumulator.cuh"
#include "lbm/boundary/dirichlet/dirichlet_eval_density.cuh"
#include "lbm/boundary/dirichlet/dirichlet_incoming_evaluation.cuh"
#include "lbm/boundary/dirichlet/dirichlet_outgoing_evaluation.cuh"

__device__ __forceinline__ void apply_boundary(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    NodeMoments &M)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    DirichletAccumulator acc{};

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        if (dir_valid(incoming_mask, i))
            dirichlet_incoming_evaluation(acc, M, pop, i);

        if (dir_valid(outgoing_mask, i))
            dirichlet_outgoing_evaluation(acc, M, i);
    }

    acc.in.normalize();

    M.mxy = (acc.rho.constant * acc.in.mxy - acc.mxy.constant) / (acc.mxy.mxy - acc.rho.mxy * acc.in.mxy);

    M.mxx = M.ux * M.ux;
    M.myy = M.uy * M.uy;

    M.rho = eval_dirichlet_density(acc, M);
}