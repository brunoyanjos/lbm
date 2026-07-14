#pragma once

#include "core/types.cuh"

namespace boundary::dirichlet
{
    struct IncomingMoments
    {
        real_t rho = r::zero;

        real_t mxx = r::zero;
        real_t mxy = r::zero;
        real_t myy = r::zero;

        __device__ __forceinline__ void normalize()
        {
            const real_t inv_rho = r::one / rho;

            mxx *= inv_rho;
            mxy *= inv_rho;
            myy *= inv_rho;
        }
    };

    struct MomentEquation
    {
        real_t constant = r::zero;

        real_t mxx = r::zero;
        real_t mxy = r::zero;
        real_t myy = r::zero;
    };

    struct Accumulator
    {
        IncomingMoments in;

        MomentEquation rho;
        MomentEquation mxx;
        MomentEquation mxy;
        MomentEquation myy;
    };
}