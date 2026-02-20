#pragma once

#include "../../core/types.cuh"
#include "../../core/indexing.cuh"
#include "../stencil_active.cuh"

__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ pop,
                                                         const LBMState &S,
                                                         int c,
                                                         int x, int y)
{
#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const int icx = Stencil::cx(i);
        const int icy = Stencil::cy(i);

        const size_t n_idx = idxGlobalPeriodic(x - icx, y - icy);

        const real_t rho = S.d_rho[c][n_idx] + RHO_0;
        const real_t ux = S.d_ux[c][n_idx];
        const real_t uy = S.d_uy[c][n_idx];
        const real_t mxx = S.d_mxx[c][n_idx];
        const real_t mxy = S.d_mxy[c][n_idx];
        const real_t myy = S.d_myy[c][n_idx];

        const real_t cx = static_cast<real_t>(icx);
        const real_t cy = static_cast<real_t>(icy);

        const real_t Hxx = cx * cx - Stencil::cs2;
        const real_t Hxy = cx * cy;
        const real_t Hyy = cy * cy - Stencil::cs2;

        pop[i] = Stencil::w(i) * rho *
                 (real_t(1) + ux * cx + uy * cy + mxx * Hxx + myy * Hyy + mxy * Hxy);
    }
}
