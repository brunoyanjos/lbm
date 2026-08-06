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

        const real_t ux_col = ux + r::half * GRAVITY_X;
        const real_t uy_col = uy + r::half * GRAVITY_Y;

        const real_t mxx_col = (r::one - OMEGA) * mxx + OMEGA * ux * ux + r::two * (r::one - OMEGA * r::half) * GRAVITY_X * ux;
        const real_t mxy_col = (r::one - OMEGA) * mxy + OMEGA * ux * uy + (r::one - OMEGA * r::half) * (GRAVITY_Y * ux + GRAVITY_X * uy);
        const real_t myy_col = (r::one - OMEGA) * myy + OMEGA * uy * uy + r::two * (r::one - OMEGA * r::half) * GRAVITY_Y * uy;

        // const real_t ux_barra = ux - r::half * GRAVITY_X;
        // const real_t uy_barra = uy - r::half * GRAVITY_Y;

        // const real_t mxx_barra = (r::one - OMEGA) * mxx + OMEGA * ux * ux + r::two * (r::one - r::half * OMEGA) * GRAVITY_X * ux;
        // const real_t mxy_barra = (r::one - OMEGA) * mxy + OMEGA * ux * uy + (r::one - r::half * OMEGA) * r::half * (GRAVITY_X * uy + GRAVITY_Y * ux);
        // const real_t myy_barra = (r::one - OMEGA) * myy + OMEGA * uy * uy + r::two * (r::one - r::half * OMEGA) * GRAVITY_Y * uy;

        pop[i] = Stencil::w(i) * rho *
                 (r::one +
                  scale_factor<MomentId::ux>() * ux_col * hermite<MomentId::ux>(i) + scale_factor<MomentId::uy>() * uy_col * hermite<MomentId::uy>(i) +
                  scale_factor<MomentId::mxx>() * mxx_col * hermite<MomentId::mxx>(i) + scale_factor<MomentId::mxy>() * mxy_col * hermite<MomentId::mxy>(i) +
                  scale_factor<MomentId::myy>() * myy_col * hermite<MomentId::myy>(i));
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
