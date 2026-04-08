#pragma once

#include "../../core/types.cuh"
#include "../../core/indexing.cuh"
#include "../stencil_active.cuh"
#include "../hermite/hermite.cuh"
#include "../moment/moment_id.cuh"

__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ pop,
                                                         const LBMState &S,
                                                         int c,
                                                         int x, int y)
{
    const real_t one_minus_omega = r::one - OMEGA;
    const real_t half_omega = r::half * OMEGA;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const int cx = Stencil::cx(i);
        const int cy = Stencil::cy(i);

        const size_t n_idx = idxGlobalPeriodic(x - cx, y - cy);

        const real_t rho = S.d_rho[c][n_idx] + RHO_0;
        const real_t ux = S.d_ux[c][n_idx];
        const real_t uy = S.d_uy[c][n_idx];
        real_t mxx = S.d_mxx[c][n_idx];
        real_t mxy = S.d_mxy[c][n_idx];
        real_t myy = S.d_myy[c][n_idx];

        mxx = one_minus_omega * mxx + half_omega * ux * ux;
        mxy = one_minus_omega * mxy + OMEGA * ux * uy;
        myy = one_minus_omega * myy + half_omega * uy * uy;

        pop[i] = Stencil::w(i) * rho *
                 (r::one +
                  ux * hermite<MomentId::ux>(i) + uy * hermite<MomentId::uy>(i) +
                  mxx * hermite<MomentId::mxx>(i) + mxy * hermite<MomentId::mxy>(i) +
                  myy * hermite<MomentId::myy>(i));
    }
}
