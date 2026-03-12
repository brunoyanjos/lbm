#pragma once

#include "../../core/types.cuh"
#include "../../core/lbm_features.cuh"
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
        const auto B = Stencil::basis(i);

        const int icx = static_cast<int>(B.cx);
        const int icy = static_cast<int>(B.cy);

        const size_t n_idx = idxGlobalPeriodic(x - icx, y - icy);

        const real_t rho = S.d_rho[c][n_idx] + Geometry::RHO_0;
        const real_t ux = S.d_ux[c][n_idx];
        const real_t uy = S.d_uy[c][n_idx];

        const real_t mxx = S.d2.mxx[c][n_idx];
        const real_t mxy = S.d2.mxy[c][n_idx];
        const real_t myy = S.d2.myy[c][n_idx];

#if LBM_HAS_REG3_CROSS
        const real_t mxxy = S.d3c.mxxy[c][n_idx];
        const real_t mxyy = S.d3c.mxyy[c][n_idx];
#endif

#if LBM_HAS_REG3_AXIAL
        const real_t mxxx = S.d3a.mxxx[c][n_idx];
        const real_t myyy = S.d3a.myyy[c][n_idx];
#endif

        pop[i] = Stencil::w(i) * rho *
                 (r::one + ux * B.cx + uy * B.cy + mxx * B.Hxx + mxy * B.Hxy + myy * B.Hyy
#if LBM_HAS_REG3_CROSS
                  + mxxy * B.Hxxy + mxyy * B.Hxyy
#endif
#if LBM_HAS_REG3_AXIAL
                  + mxxx * B.Hxxx + myyy * B.Hyyy
#endif
                 );
    }
}
