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

__device__ __forceinline__ real_t rel_step(real_t x_new, real_t x_old)
{
    const real_t eps = real_t(1e-12);
    const real_t denom = fmax(r_abs(x_new), eps);
    return r_abs(x_new - x_old) / denom;
}

__device__ __forceinline__ real_t eval_row(
    const real_t *__restrict__ A,
    real_t B,
    real_t ux, real_t uy,
    real_t mxx, real_t mxy, real_t myy)
{
    const real_t ux2 = ux * ux;
    const real_t uxuy = ux * uy;
    const real_t uy2 = uy * uy;

    return ux * A[0] + uy * A[1] +
           ux2 * A[2] + uxuy * A[3] + uy2 * A[4] +
           mxx * A[5] + mxy * A[6] + myy * A[7] -
           B;
}

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

    build_fluid_system(S, M.ux, M.uy, acc);

    SystemData<UnknownMomentList<true>::size> G{};

    real_t error = r::one;
    int it = r::zero;
    const int it_max = 50;

    while (error > real_t(1e-6) && it++ < it_max)
    {
        const real_t ux_old = M.ux;
        const real_t uy_old = M.uy;
        const real_t mxx_old = M.mxx;
        const real_t mxy_old = M.mxy;
        const real_t myy_old = M.myy;

        for (int i = 0; i < 5; ++i)
        {
            G.coeff(i, 0) = S.coeff(i, 0) + r::two * S.coeff(i, 2) * M.ux + S.coeff(i, 3) * M.uy;
            G.coeff(i, 1) = S.coeff(i, 1) + S.coeff(i, 3) * M.ux + r::two * S.coeff(i, 4) * M.uy;
            G.coeff(i, 2) = S.coeff(i, 5);
            G.coeff(i, 3) = S.coeff(i, 6);
            G.coeff(i, 4) = S.coeff(i, 7);

            G.b[i] = -eval_row(&S.A[i * S.cols], S.b[i], M.ux, M.uy, M.mxx, M.mxy, M.myy);
        }

        gaussianElimination(G);

        M.ux = ux_old + G.x[0];
        M.uy = uy_old + G.x[1];
        M.mxx = mxx_old + G.x[2];
        M.mxy = mxy_old + G.x[3];
        M.myy = myy_old + G.x[4];

        error = rel_step(M.ux, ux_old);
        error = fmax(error, rel_step(M.uy, uy_old));
        error = fmax(error, rel_step(M.mxx, mxx_old));
        error = fmax(error, rel_step(M.mxy, mxy_old));
        error = fmax(error, rel_step(M.myy, myy_old));
    }

    M.rho = eval_density(acc, M);
}
