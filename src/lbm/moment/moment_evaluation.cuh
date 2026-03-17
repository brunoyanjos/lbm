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
        const real_t cx = static_cast<real_t>(Stencil::cx(i));
        const real_t cy = static_cast<real_t>(Stencil::cy(i));

        const real_t Hxx = cx * cx - Stencil::cs2;
        const real_t Hxy = cx * cy;
        const real_t Hyy = cy * cy - Stencil::cs2;

        const real_t fi = pop[i];

        M.rho += fi;

        M.ux += fi * cx;
        M.uy += fi * cy;

        M.mxx += fi * Hxx;
        M.mxy += fi * Hxy;
        M.myy += fi * Hyy;
    }

    const real_t inv_rho = real_t(1) / M.rho;

    M.ux *= inv_rho;
    M.uy *= inv_rho;
    M.mxx *= inv_rho;
    M.mxy *= inv_rho;
    M.myy *= inv_rho;
}