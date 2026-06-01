#pragma once

#include "../../core/types.cuh"
#include "../../core/indexing.cuh"
#include "../../core/simulation_config.h"
#include "../../core/physics.h"
#include "../stencil_active.cuh"
#include "../hermite/hermite.cuh"
#include "../moment/moment_id.cuh"
#include "lbm/state/lbm_state.cuh"

template <int RegOrder, bool Rec, bool HighOrder>
__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ pop,
                                                         const LBMStateFor<RegOrder, Rec, HighOrder> &S,
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

template <bool HighOrder>
__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ pop,
                                                         const LBMStateFor<3, true, HighOrder> &S,
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
        const real_t mxxy = ux * mxy + uy * mxx - ux * ux * uy;
        const real_t mxyy = uy * mxy + ux * myy - ux * uy * uy;

        real_t expansion =
            r::one +
            ux * hermite<MomentId::ux>(i) + uy * hermite<MomentId::uy>(i) +
            mxx * hermite<MomentId::mxx>(i) + mxy * hermite<MomentId::mxy>(i) +
            myy * hermite<MomentId::myy>(i) + mxxy * hermite<MomentId::mxxy>(i) +
            mxyy * hermite<MomentId::mxyy>(i);

        if constexpr (HighOrder)
        {
            const real_t mxxx = ux * mxx - r::third * ux * ux * ux;
            const real_t myyy = uy * myy - r::third * uy * uy * uy;
            expansion += mxxx * hermite<MomentId::mxxx>(i) + myyy * hermite<MomentId::myyy>(i);
        }

        pop[i] = Stencil::w(i) * rho * expansion;
    }
}

template <bool HighOrder>
__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ pop,
                                                         const LBMStateFor<3, false, HighOrder> &S,
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
        const real_t mxxy = S.d_mxxy[c][n_idx];
        const real_t mxyy = S.d_mxyy[c][n_idx];

        real_t expansion =
            r::one +
            ux * hermite<MomentId::ux>(i) + uy * hermite<MomentId::uy>(i) +
            mxx * hermite<MomentId::mxx>(i) + mxy * hermite<MomentId::mxy>(i) +
            myy * hermite<MomentId::myy>(i) + mxxy * hermite<MomentId::mxxy>(i) +
            mxyy * hermite<MomentId::mxyy>(i);

        if constexpr (HighOrder)
        {
            const real_t mxxx = S.d_mxxx[c][n_idx];
            const real_t myyy = S.d_myyy[c][n_idx];
            expansion += mxxx * hermite<MomentId::mxxx>(i) + myyy * hermite<MomentId::myyy>(i);
        }

        pop[i] = Stencil::w(i) * rho * expansion;
    }
}
