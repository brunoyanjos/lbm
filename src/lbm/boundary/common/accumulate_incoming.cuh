#pragma once

#include "core/types.cuh"

#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/moment_id.cuh"

#include "lbm/boundary/common/accumulator.cuh"
#include "lbm/boundary/common/factor.cuh"

template <bool HasVelocity>
__device__ __forceinline__ void accumulate_incoming(MomentAccumulator<HasVelocity> &acc, const real_t *__restrict__ pop, int i)
{
    acc.in.rho += pop[i];

    if constexpr (HasVelocity)
    {
        acc.in.ux += pop[i] * hermite<MomentId::ux>(i);
        acc.in.uy += pop[i] * hermite<MomentId::uy>(i);
    }

    acc.in.mxx += pop[i] * hermite<MomentId::mxx>(i);
    acc.in.mxy += pop[i] * hermite<MomentId::mxy>(i);
    acc.in.myy += pop[i] * hermite<MomentId::myy>(i);

    if constexpr (HasVelocity)
    {
        acc.ux.rho += moment_factor<MomentId::rho, MomentId::ux>(i);
        acc.ux.ux += moment_factor<MomentId::ux, MomentId::ux>(i);
        acc.ux.uy += moment_factor<MomentId::uy, MomentId::ux>(i);
        acc.ux.mxx += moment_factor<MomentId::mxx, MomentId::ux>(i);
        acc.ux.mxy += moment_factor<MomentId::mxy, MomentId::ux>(i);
        acc.ux.myy += moment_factor<MomentId::myy, MomentId::ux>(i);

        acc.uy.rho += moment_factor<MomentId::rho, MomentId::uy>(i);
        acc.uy.ux += moment_factor<MomentId::ux, MomentId::uy>(i);
        acc.uy.uy += moment_factor<MomentId::uy, MomentId::uy>(i);
        acc.uy.mxx += moment_factor<MomentId::mxx, MomentId::uy>(i);
        acc.uy.mxy += moment_factor<MomentId::mxy, MomentId::uy>(i);
        acc.uy.myy += moment_factor<MomentId::myy, MomentId::uy>(i);
    }

    acc.mxx.rho += moment_factor<MomentId::rho, MomentId::mxx>(i);
    acc.mxx.ux += moment_factor<MomentId::ux, MomentId::mxx>(i);
    acc.mxx.uy += moment_factor<MomentId::uy, MomentId::mxx>(i);
    acc.mxx.mxx += moment_factor<MomentId::mxx, MomentId::mxx>(i);
    acc.mxx.mxy += moment_factor<MomentId::mxy, MomentId::mxx>(i);
    acc.mxx.myy += moment_factor<MomentId::myy, MomentId::mxx>(i);

    acc.mxy.rho += moment_factor<MomentId::rho, MomentId::mxy>(i);
    acc.mxy.ux += moment_factor<MomentId::ux, MomentId::mxy>(i);
    acc.mxy.uy += moment_factor<MomentId::uy, MomentId::mxy>(i);
    acc.mxy.mxx += moment_factor<MomentId::mxx, MomentId::mxy>(i);
    acc.mxy.mxy += moment_factor<MomentId::mxy, MomentId::mxy>(i);
    acc.mxy.myy += moment_factor<MomentId::myy, MomentId::mxy>(i);

    acc.myy.rho += moment_factor<MomentId::rho, MomentId::myy>(i);
    acc.myy.ux += moment_factor<MomentId::ux, MomentId::myy>(i);
    acc.myy.uy += moment_factor<MomentId::uy, MomentId::myy>(i);
    acc.myy.mxx += moment_factor<MomentId::mxx, MomentId::myy>(i);
    acc.myy.mxy += moment_factor<MomentId::mxy, MomentId::myy>(i);
    acc.myy.myy += moment_factor<MomentId::myy, MomentId::myy>(i);
}