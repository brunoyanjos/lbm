#pragma once

#include "core/types.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/moment_id.cuh"
#include "lbm/moment/node_moments.cuh"

template <int Order, bool Rec, bool HighOrder>
__device__ __forceinline__ void evaluate_moments_from_pop(const real_t *__restrict__ pop,
                                                          NodeMomentsFor<Order, Rec, HighOrder> &M)
{
    M.rho = r::zero;
    M.ux = r::zero;
    M.uy = r::zero;
    M.mxx = r::zero;
    M.mxy = r::zero;
    M.myy = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const real_t fi = (pop[i] + Stencil::w(i));

        M.rho += fi;

        M.ux += fi * hermite<MomentId::ux>(i);
        M.uy += fi * hermite<MomentId::uy>(i);

        M.mxx += fi * hermite<MomentId::mxx>(i);
        M.mxy += fi * hermite<MomentId::mxy>(i);
        M.myy += fi * hermite<MomentId::myy>(i);
    }

    const real_t inv_rho = r::one / M.rho;

    M.ux *= inv_rho;
    M.uy *= inv_rho;
    M.mxx *= inv_rho;
    M.mxy *= inv_rho;
    M.myy *= inv_rho;
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
        const real_t fi = pop[i] + Stencil::w(i);

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
