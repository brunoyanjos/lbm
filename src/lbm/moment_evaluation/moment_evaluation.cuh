#pragma once

#include "../../core/types.cuh"
#include "../stencil_active.cuh"

__device__ __forceinline__ void evaluate_moments_from_pop(const real_t *__restrict__ pop,
                                                          real_t &rho,
                                                          real_t &ux, real_t &uy,
                                                          real_t &mxx, real_t &mxy, real_t &myy)
{
    rho = r::zero;
    ux = r::zero;
    uy = r::zero;
    mxx = r::zero;
    mxy = r::zero;
    myy = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const auto B = Stencil::basis(i);

        const real_t fi = pop[i];

        rho += fi;

        ux += fi * B.cx;
        uy += fi * B.cy;

        mxx += fi * B.Hxx;
        mxy += fi * B.Hxy;
        myy += fi * B.Hyy;
    }

    const real_t inv_rho = r::one / rho;

    ux *= inv_rho;
    uy *= inv_rho;
    mxx *= inv_rho;
    mxy *= inv_rho;
    myy *= inv_rho;
}