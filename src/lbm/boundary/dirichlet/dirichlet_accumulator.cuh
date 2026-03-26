#pragma once

#include "core/types.cuh"

struct DirichletIncomingMoments
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

struct DirichletMomentEquation
{
    real_t constant = r::zero;

    real_t mxx = r::zero;
    real_t mxy = r::zero;
    real_t myy = r::zero;
};

struct DirichletAccumulator
{
    DirichletIncomingMoments in;

    DirichletMomentEquation rho;

    DirichletMomentEquation mxx;
    DirichletMomentEquation mxy;
    DirichletMomentEquation myy;
};