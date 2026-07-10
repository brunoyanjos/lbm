#pragma once

#include "core/types.cuh"

namespace boundary::fluid
{
    struct IncomingMomentAccumulator
    {
        real_t rho = r::zero;

        real_t ux = r::zero;
        real_t uy = r::zero;

        real_t mxx = r::zero;
        real_t mxy = r::zero;
        real_t myy = r::zero;

        __device__ __forceinline__ void normalize()
        {
            const real_t inv_rho = r::one / rho;

            ux *= inv_rho;
            uy *= inv_rho;
            mxx *= inv_rho;
            mxy *= inv_rho;
            myy *= inv_rho;
        }
    };

    struct DensityEqAccumulator
    {
        real_t rho = r::zero;

        real_t ux = r::zero;
        real_t uy = r::zero;

        real_t uxux = r::zero;
        real_t uxuy = r::zero;
        real_t uyuy = r::zero;

        real_t mxx = r::zero;
        real_t mxy = r::zero;
        real_t myy = r::zero;
    };

    struct SecondOrderMomentEqAccumulator
    {
        real_t rho = r::zero;

        real_t ux = r::zero;
        real_t uy = r::zero;

        real_t mxx = r::zero;
        real_t mxy = r::zero;
        real_t myy = r::zero;
    };

    struct Accumulator
    {
        IncomingMomentAccumulator in;

        DensityEqAccumulator rho;

        SecondOrderMomentEqAccumulator ux;
        SecondOrderMomentEqAccumulator uy;

        SecondOrderMomentEqAccumulator mxx;
        SecondOrderMomentEqAccumulator mxy;
        SecondOrderMomentEqAccumulator myy;
    };
}