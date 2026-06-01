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

    struct DensityEqAccumulatorSecondOrder
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

    template <int RegOrder, bool Rec>
    struct DensityEqAccumulator : DensityEqAccumulatorSecondOrder
    {
    };

    template <>
    struct DensityEqAccumulator<3, false> : DensityEqAccumulatorSecondOrder
    {
        real_t uxuxux = r::zero;
        real_t uxuxuy = r::zero;
        real_t uxuyuy = r::zero;
        real_t uyuyuy = r::zero;
    };

    template <>
    struct DensityEqAccumulator<3, true> : DensityEqAccumulator<3, false>
    {
        real_t uxmxx = r::zero;
        real_t uymxx = r::zero;
        real_t uxmyy = r::zero;
        real_t uymyy = r::zero;
        real_t uxmxy = r::zero;
        real_t uymxy = r::zero;
    };

    struct SecondOrderMomentEqAccumulatorSecondOrder
    {
        real_t rho = r::zero;

        real_t ux = r::zero;
        real_t uy = r::zero;

        real_t mxx = r::zero;
        real_t mxy = r::zero;
        real_t myy = r::zero;
    };

    template <int RegOrder, bool Rec>
    struct SecondOrderMomentEqAccumulator : SecondOrderMomentEqAccumulatorSecondOrder
    {
    };

    template <>
    struct SecondOrderMomentEqAccumulator<3, false> : SecondOrderMomentEqAccumulatorSecondOrder
    {
        real_t uxuxux = r::zero;
        real_t uxuxuy = r::zero;
        real_t uxuyuy = r::zero;
        real_t uyuyuy = r::zero;
    };

    template <>
    struct SecondOrderMomentEqAccumulator<3, true> : SecondOrderMomentEqAccumulator<3, false>
    {
        real_t uxmxx = r::zero;
        real_t uymxx = r::zero;
        real_t uxmyy = r::zero;
        real_t uymyy = r::zero;
        real_t uxmxy = r::zero;
        real_t uymxy = r::zero;
    };

    template <int RegOrder, bool Rec>
    struct Accumulator
    {
        IncomingMomentAccumulator in;

        DensityEqAccumulator<RegOrder, Rec> rho;

        SecondOrderMomentEqAccumulator<RegOrder, Rec> ux;
        SecondOrderMomentEqAccumulator<RegOrder, Rec> uy;

        SecondOrderMomentEqAccumulator<RegOrder, Rec> mxx;
        SecondOrderMomentEqAccumulator<RegOrder, Rec> mxy;
        SecondOrderMomentEqAccumulator<RegOrder, Rec> myy;
    };
}