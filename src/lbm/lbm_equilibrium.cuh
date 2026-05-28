#pragma once

#include "../core/types.cuh"
#include "../core/simulation_config.h"
#include "hermite/hermite.cuh"
#include "moment/moment_id.cuh"
#include "moment/scale_factor.cuh"
#include "stencil_active.cuh"

template <int RegOrder, bool HighOrder>
__device__ __forceinline__ real_t equilibrium_expansion(int i, real_t ux, real_t uy)
{
    real_t e =
        r::one +
        scale_factor<MomentId::ux>() * ux * hermite<MomentId::ux>(i) +
        scale_factor<MomentId::uy>() * uy * hermite<MomentId::uy>(i) +
        scale_factor<MomentId::mxx>() * ux * ux * hermite<MomentId::mxx>(i) +
        scale_factor<MomentId::mxy>() * ux * uy * hermite<MomentId::mxy>(i) +
        scale_factor<MomentId::myy>() * uy * uy * hermite<MomentId::myy>(i);

    if constexpr (RegOrder >= 3)
    {
        e += scale_factor<MomentId::mxxy>() * ux * ux * uy * hermite<MomentId::mxxy>(i);
        e += scale_factor<MomentId::mxyy>() * ux * uy * uy * hermite<MomentId::mxyy>(i);

        if constexpr (HighOrder)
        {
            e += scale_factor<MomentId::mxxx>() * ux * ux * ux * hermite<MomentId::mxxx>(i);
            e += scale_factor<MomentId::myyy>() * uy * uy * uy * hermite<MomentId::myyy>(i);
        }
    }

    return e;
}

__device__ __forceinline__ void equilibrium(real_t *pop, real_t rho, real_t ux, real_t uy)
{
#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
        pop[i] = Stencil::w(i) * rho *
                 equilibrium_expansion<REG_ORDER, Stencil::high_order>(i, ux, uy);
}
