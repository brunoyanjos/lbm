#pragma once

#include "../../core/types.cuh"
#include "../../core/indexing.cuh"
#include "../../core/simulation_config.h"
#include "../../core/physics.h"
#include "../stencil_active.cuh"
#include "../hermite/hermite.cuh"
#include "../moment/moment_id.cuh"
#include "lbm/state/lbm_state.cuh"

__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ pop,
                                                         const LBMState &S,
                                                         int c, int x, int y)
{
#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const int cx = Stencil::cx(i);
        const int cy = Stencil::cy(i);

        const size_t n_idx = idxGlobalPeriodic(x - cx, y - cy);

        const real_t rho = S.d_rho[c][n_idx] + RHO_0;
        const real_t ux = S.d_ux[c][n_idx];
        const real_t uy = S.d_uy[c][n_idx];
        const real_t mxx = S.d_mxx[c][n_idx];
        const real_t mxy = S.d_mxy[c][n_idx];
        const real_t myy = S.d_myy[c][n_idx];

        pop[i] = Stencil::w(i) * rho *
                 (r::one +
                  ux * hermite<MomentId::ux>(i) + uy * hermite<MomentId::uy>(i) +
                  mxx * hermite<MomentId::mxx>(i) + mxy * hermite<MomentId::mxy>(i) +
                  myy * hermite<MomentId::myy>(i));
    }
}
