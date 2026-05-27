#pragma once

#include "core/types.cuh"

template <bool HasVelocity>
struct IncomingMomentAccumulator;

template <>
struct IncomingMomentAccumulator<true>
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

template <>
struct IncomingMomentAccumulator<false>
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

struct LinearMomentEqAccumulator
{
    real_t rho = r::zero;

    real_t ux = r::zero;
    real_t uy = r::zero;

    real_t mxx = r::zero;
    real_t mxy = r::zero;
    real_t myy = r::zero;
};

template <bool HasVelocity>
struct MomentAccumulator;

template <>
struct MomentAccumulator<true>
{
    IncomingMomentAccumulator<true> in;
    DensityEqAccumulator rho;
    LinearMomentEqAccumulator ux;
    LinearMomentEqAccumulator uy;
    LinearMomentEqAccumulator mxx;
    LinearMomentEqAccumulator mxy;
    LinearMomentEqAccumulator myy;
};

template <>
struct MomentAccumulator<false>
{
    IncomingMomentAccumulator<false> in;
    DensityEqAccumulator rho;
    LinearMomentEqAccumulator mxx;
    LinearMomentEqAccumulator mxy;
    LinearMomentEqAccumulator myy;
};