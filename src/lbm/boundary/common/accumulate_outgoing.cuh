#pragma once

#include "core/physics.h"
#include "core/types.cuh"

#include "lbm/moment/moment_id.cuh"

#include "lbm/boundary/common/accumulator.cuh"
#include "lbm/boundary/common/factor.cuh"

template <bool HasVelocity>
__device__ __forceinline__ void accumulate_outgoing(MomentAccumulator<HasVelocity> &acc, int i)
{
    acc.rho.rho += moment_factor<MomentId::rho, MomentId::rho>(i);
    acc.rho.ux += moment_factor<MomentId::ux, MomentId::rho>(i);
    acc.rho.uy += moment_factor<MomentId::uy, MomentId::rho>(i);
    acc.rho.uxux += OMEGA * moment_factor<MomentId::mxx, MomentId::rho>(i);
    acc.rho.uxuy += OMEGA * moment_factor<MomentId::mxy, MomentId::rho>(i);
    acc.rho.uyuy += OMEGA * moment_factor<MomentId::myy, MomentId::rho>(i);
    acc.rho.mxx += (r::one - OMEGA) * moment_factor<MomentId::mxx, MomentId::rho>(i);
    acc.rho.mxy += (r::one - OMEGA) * moment_factor<MomentId::mxy, MomentId::rho>(i);
    acc.rho.myy += (r::one - OMEGA) * moment_factor<MomentId::myy, MomentId::rho>(i);
}