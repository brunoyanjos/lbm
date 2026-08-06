#pragma once

#include "../../core/types.cuh"
#include "../../core/indexing.cuh"
#include "../../core/simulation_config.h"
#include "../../core/physics.h"
#include "../stencil_active.cuh"
#include "../hermite/hermite.cuh"
#include "../moment/moment_id.cuh"
#include "lbm/state/lbm_state.cuh"

struct Normal
{
    real_t x = r::zero;
    real_t y = r::zero;
};

template <int RegOrder, bool Rec, bool HighOrder>
[[nodiscard]]
__device__ __forceinline__ Normal evaluate_normal(int c, int x, int y,
                                                  const LBMStateFor<RegOrder, Rec, HighOrder> &S)
{
    Normal n{};

    n.x = r::zero;
    n.y = r::zero;

    for (int i = 0; i < Stencil::Q; ++i)
    {
        const int cx = Stencil::cx(i);
        const int cy = Stencil::cy(i);

        const size_t n_idx = idxGlobalPeriodic(x + cx, y + cy);

        const real_t rhoA = S.d_rhoA[c][n_idx];
        const real_t rhoB = S.d_rhoB[c][n_idx];

        const real_t phi = (rhoA - rhoB) / (rhoA + rhoB);

        n.x += Stencil::w(i) * phi * cx * Stencil::as2;
        n.y += Stencil::w(i) * phi * cy * Stencil::as2;
    }

    const real_t norm = dsqrt(n.x * n.x + n.y * n.y);

    const real_t inv_norm = norm > r::zero ? r::one / norm : r::zero;

    n.x *= inv_norm;
    n.y *= inv_norm;

    return n;
}

template <int RegOrder, bool Rec, bool HighOrder>
__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ popA, real_t *__restrict__ popB,
                                                         const LBMStateFor<RegOrder, Rec, HighOrder> &S,
                                                         int c, int x, int y)
{
#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const int cx = Stencil::cx(i);
        const int cy = Stencil::cy(i);

        const size_t n_idx = idxGlobalPeriodic(x - cx, y - cy);

        const real_t rhoA = S.d_rhoA[c][n_idx];
        const real_t uxA = S.d_uxA[c][n_idx];
        const real_t uyA = S.d_uyA[c][n_idx];
        const real_t mxxA = S.d_mxxA[c][n_idx];
        const real_t mxyA = S.d_mxyA[c][n_idx];
        const real_t myyA = S.d_myyA[c][n_idx];

        const real_t rhoB = S.d_rhoB[c][n_idx];
        const real_t uxB = S.d_uxB[c][n_idx];
        const real_t uyB = S.d_uyB[c][n_idx];
        const real_t mxxB = S.d_mxxB[c][n_idx];
        const real_t mxyB = S.d_mxyB[c][n_idx];
        const real_t myyB = S.d_myyB[c][n_idx];

        const real_t rho = rhoA + rhoB;

        const real_t inv_rho = r::one / rho;
        const real_t xA = rhoA * inv_rho;
        const real_t xB = rhoB * inv_rho;

        const real_t TAU_A = TAU;
        const real_t TAU_B = TAU * 1000000;

        const real_t OMEGA_A = static_cast<real_t>(1) / TAU_A;
        const real_t OMEGA_B = static_cast<real_t>(1) / TAU_B;

        const real_t OMEGA_MIX = xA * OMEGA_A + xB * OMEGA_B;

        const Normal n = evaluate_normal(c, x - cx, y - cy, S);

        const real_t uxA_col = uxA + r::half * GRAVITY_X + BETA * xB * n.x * Stencil::cs2;
        const real_t uyA_col = uyA + r::half * GRAVITY_Y + BETA * xB * n.y * Stencil::cs2;

        const real_t mxxA_col = (r::one - OMEGA_MIX) * mxxA + OMEGA_MIX * uxA * uxA + r::two * (r::one - OMEGA_MIX * r::half) * GRAVITY_X * uxA;
        const real_t mxyA_col = (r::one - OMEGA_MIX) * mxyA + OMEGA_MIX * uxA * uyA + (r::one - OMEGA_MIX * r::half) * (GRAVITY_Y * uxA + GRAVITY_X * uyA);
        const real_t myyA_col = (r::one - OMEGA_MIX) * myyA + OMEGA_MIX * uyA * uyA + r::two * (r::one - OMEGA_MIX * r::half) * GRAVITY_Y * uyA;

        const real_t uxB_col = uxB + r::half * GRAVITY_X - BETA * xA * n.x * Stencil::cs2;
        const real_t uyB_col = uyB + r::half * GRAVITY_Y - BETA * xA * n.y * Stencil::cs2;

        const real_t mxxB_col = (r::one - OMEGA_MIX) * mxxB + OMEGA_MIX * uxB * uxB + r::two * (r::one - OMEGA_MIX * r::half) * GRAVITY_X * uxB;
        const real_t mxyB_col = (r::one - OMEGA_MIX) * mxyB + OMEGA_MIX * uxB * uyB + (r::one - OMEGA_MIX * r::half) * (GRAVITY_Y * uxB + GRAVITY_X * uyB);
        const real_t myyB_col = (r::one - OMEGA_MIX) * myyB + OMEGA_MIX * uyB * uyB + r::two * (r::one - OMEGA_MIX * r::half) * GRAVITY_Y * uyB;

        popA[i] = Stencil::w(i) * rhoA *
                  (r::one +
                   scale_factor<MomentId::ux>() * uxA_col * hermite<MomentId::ux>(i) + scale_factor<MomentId::uy>() * uyA_col * hermite<MomentId::uy>(i) +
                   scale_factor<MomentId::mxx>() * mxxA_col * hermite<MomentId::mxx>(i) + scale_factor<MomentId::mxy>() * mxyA_col * hermite<MomentId::mxy>(i) +
                   scale_factor<MomentId::myy>() * myyA_col * hermite<MomentId::myy>(i));

        popB[i] = Stencil::w(i) * rhoB *
                  (r::one +
                   scale_factor<MomentId::ux>() * uxB_col * hermite<MomentId::ux>(i) + scale_factor<MomentId::uy>() * uyB_col * hermite<MomentId::uy>(i) +
                   scale_factor<MomentId::mxx>() * mxxB_col * hermite<MomentId::mxx>(i) + scale_factor<MomentId::mxy>() * mxyB_col * hermite<MomentId::mxy>(i) +
                   scale_factor<MomentId::myy>() * myyB_col * hermite<MomentId::myy>(i));
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

        const real_t rho = S.d_rho[c][n_idx] + RHOA_0;
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

        const real_t rho = S.d_rho[c][n_idx] + RHOA_0;
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
