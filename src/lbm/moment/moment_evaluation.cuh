#pragma once

#include "core/types.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/moment_id.cuh"
#include "lbm/moment/node_moments.cuh"

template <int Order, bool Rec, bool HighOrder>
__device__ __forceinline__ void evaluate_moments_from_pop(const real_t *__restrict__ popA,
                                                          const real_t *__restrict__ popB,
                                                          NodeMomentsFor<Order, Rec, HighOrder> &M)
{
    M.rhoA = r::zero;
    M.uxA = r::zero;
    M.uyA = r::zero;
    M.mxxA = r::zero;
    M.mxyA = r::zero;
    M.myyA = r::zero;

    M.rhoB = r::zero;
    M.uxB = r::zero;
    M.uyB = r::zero;
    M.mxxB = r::zero;
    M.mxyB = r::zero;
    M.myyB = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const real_t fiA = popA[i];
        const real_t fiB = popB[i];

        M.rhoA += fiA;
        M.rhoB += fiB;

        M.uxA += (fiA + fiB) * hermite<MomentId::ux>(i);
        M.uyA += (fiA + fiB) * hermite<MomentId::uy>(i);

        M.uxB += (fiA + fiB) * hermite<MomentId::ux>(i);
        M.uyB += (fiA + fiB) * hermite<MomentId::uy>(i);

        M.mxxA += (fiA + fiB) * hermite<MomentId::mxx>(i);
        M.mxyA += (fiA + fiB) * hermite<MomentId::mxy>(i);
        M.myyA += (fiA + fiB) * hermite<MomentId::myy>(i);

        M.mxxB += (fiA + fiB) * hermite<MomentId::mxx>(i);
        M.mxyB += (fiA + fiB) * hermite<MomentId::mxy>(i);
        M.myyB += (fiA + fiB) * hermite<MomentId::myy>(i);
    }

    const real_t inv_rho = r::one / (M.rhoA + M.rhoB);

    M.uxA *= inv_rho;
    M.uyA *= inv_rho;
    M.mxxA *= inv_rho;
    M.mxyA *= inv_rho;
    M.myyA *= inv_rho;

    M.uxB *= inv_rho;
    M.uyB *= inv_rho;
    M.mxxB *= inv_rho;
    M.mxyB *= inv_rho;
    M.myyB *= inv_rho;

    M.uxA += GRAVITY_X * r::half;
    M.uyA += GRAVITY_Y * r::half;

    M.mxxA += GRAVITY_X * M.uxA;
    M.mxyA += r::half * (GRAVITY_Y * M.uxA + GRAVITY_X * M.uyA);
    M.myyA += GRAVITY_Y * M.uyA;

    M.uxB += GRAVITY_X * r::half;
    M.uyB += GRAVITY_Y * r::half;

    M.mxxB += GRAVITY_X * M.uxB;
    M.mxyB += r::half * (GRAVITY_Y * M.uxB + GRAVITY_X * M.uyB);
    M.myyB += GRAVITY_Y * M.uyB;
}

template <bool HighOrder>
__device__ __forceinline__ void evaluate_moments_from_pop(const real_t *__restrict__ pop,
                                                          NodeMomentsFor<3, false, HighOrder> &M)
{
    M.rho = r::zero;
    M.ux = r::zero;
    M.uy = r::zero;
    M.mxx = r::zero;
    M.mxy = r::zero;
    M.myy = r::zero;
    M.mxxy = r::zero;
    M.mxyy = r::zero;
    if constexpr (HighOrder)
    {
        M.mxxx = r::zero;
        M.myyy = r::zero;
    }

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const real_t fi = pop[i];

        M.rho += fi;

        M.ux += fi * hermite<MomentId::ux>(i);
        M.uy += fi * hermite<MomentId::uy>(i);

        M.mxx += fi * hermite<MomentId::mxx>(i);
        M.mxy += fi * hermite<MomentId::mxy>(i);
        M.myy += fi * hermite<MomentId::myy>(i);

        M.mxxy += fi * hermite<MomentId::mxxy>(i);
        M.mxyy += fi * hermite<MomentId::mxyy>(i);
        if constexpr (HighOrder)
        {
            M.mxxx += fi * hermite<MomentId::mxxx>(i);
            M.myyy += fi * hermite<MomentId::myyy>(i);
        }
    }

    const real_t inv_rho = r::one / M.rho;

    M.ux *= inv_rho;
    M.uy *= inv_rho;
    M.mxx *= inv_rho;
    M.mxy *= inv_rho;
    M.myy *= inv_rho;
    M.mxxy *= inv_rho;
    M.mxyy *= inv_rho;
    if constexpr (HighOrder)
    {
        M.mxxx *= inv_rho;
        M.myyy *= inv_rho;
    }
}