#pragma once

#include "core/types.cuh"

namespace boundary::dirichlet
{
    struct IncomingMoments
    {
        real_t rho = r::zero;

        real_t mxy = r::zero;

        __device__ __forceinline__ void normalize()
        {
            const real_t inv_rho = r::one / rho;

            mxy *= inv_rho;
        }
    };

    struct MomentEquation
    {
        real_t constant = r::zero;

        real_t mxy = r::zero;
    };

    struct Accumulator
    {
        IncomingMoments in;

        MomentEquation rho;
        MomentEquation mxy;
    };
}