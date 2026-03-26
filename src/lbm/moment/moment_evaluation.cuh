#pragma once

#include "../../core/types.cuh"
#include "../stencil_active.cuh"

__device__ __forceinline__ void evaluate_moments_from_pop(const real_t *__restrict__ pop,
                                                          NodeMoments &M)
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
        const real_t fi = pop[i];

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