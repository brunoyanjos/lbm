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
    __device__ __forceinline__ void apply_boundary(
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
                evaluate_incoming(acc, M, pop, i);

            if (dir_valid(outgoing_mask, i))
                evaluate_outgoing(acc, M, i);
        }

        acc.in.normalize();

        SystemData<3> S;

        // mxx
        S.coeff(0, 0) = acc.mxx.mxx - acc.rho.mxx * acc.in.mxx; // mxx
        S.coeff(0, 1) = acc.mxx.mxy - acc.rho.mxy * acc.in.mxx; // mxy
        S.coeff(0, 2) = acc.mxx.myy - acc.rho.myy * acc.in.mxx; // myy

        S.b[0] = acc.rho.constant * acc.in.mxx - acc.mxx.constant;

        // mxx
        S.coeff(1, 0) = acc.mxy.mxx - acc.rho.mxx * acc.in.mxy; // mxx
        S.coeff(1, 1) = acc.mxy.mxy - acc.rho.mxy * acc.in.mxy; // mxy
        S.coeff(1, 2) = acc.mxy.myy - acc.rho.myy * acc.in.mxy; // myy

        S.b[1] = acc.rho.constant * acc.in.mxy - acc.mxy.constant;

        // myy
        S.coeff(2, 0) = acc.myy.mxx - acc.rho.mxx * acc.in.myy; // mxx
        S.coeff(2, 1) = acc.myy.mxy - acc.rho.mxy * acc.in.myy; // mxy
        S.coeff(2, 2) = acc.myy.myy - acc.rho.myy * acc.in.myy; // myy

        S.b[2] = acc.rho.constant * acc.in.myy - acc.myy.constant;

        gaussianElimination(S);

        M.mxx = S.x[0];
        M.mxy = S.x[1];
        M.myy = S.x[2];

        M.rho = eval_density(acc, M);
    }
}