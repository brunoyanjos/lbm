#pragma once

#include "core/types.cuh"

struct DirichletIncomingMoments
{
    real_t rho = r::zero;

    real_t mxy = r::zero;

    __device__ __forceinline__ void normalize()
    {
        const real_t inv_rho = r::one / rho;

        mxy *= inv_rho;
    }
};

struct DirichletMomentEquation
{
    real_t constant = r::zero;

    real_t mxy = r::zero;
};

struct DirichletAccumulator
{
    DirichletIncomingMoments in;

    DirichletMomentEquation rho;
    DirichletMomentEquation mxy;
};