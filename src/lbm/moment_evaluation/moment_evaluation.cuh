#pragma once

#include "../../core/types.cuh"
#include "../stencil_active.cuh"

__device__ __forceinline__ void evaluate_moments_from_pop(const real_t *__restrict__ pop,
                                                          real_t &rho,
                                                          real_t &ux, real_t &uy,
                                                          real_t &mxx, real_t &mxy, real_t &myy)
{
    rho = real_t(0);
    ux = real_t(0);
    uy = real_t(0);
    mxx = real_t(0);
    mxy = real_t(0);
    myy = real_t(0);

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const real_t cx = static_cast<real_t>(Stencil::cx(i));
        const real_t cy = static_cast<real_t>(Stencil::cy(i));

        const real_t Hxx = cx * cx - Stencil::cs2;
        const real_t Hxy = cx * cy;
        const real_t Hyy = cy * cy - Stencil::cs2;

        const real_t fi = pop[i];

        rho += fi;

        ux += fi * cx;
        uy += fi * cy;

        mxx += fi * Hxx;
        mxy += fi * Hxy;
        myy += fi * Hyy;
    }

    const real_t inv_rho = real_t(1) / rho;

    ux *= inv_rho;
    uy *= inv_rho;
    mxx *= inv_rho;
    mxy *= inv_rho;
    myy *= inv_rho;
}