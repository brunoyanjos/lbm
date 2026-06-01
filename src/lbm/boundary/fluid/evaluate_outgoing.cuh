#pragma once

#include "core/physics.h"
#include "core/types.cuh"

#include "lbm/moment/moment_id.cuh"

#include "lbm/boundary/fluid/accumulator.cuh"
#include "lbm/boundary/common/factor.cuh"

namespace boundary::fluid
{
    template <int RegOrder, bool Rec>
    __device__ __forceinline__ void evaluate_outgoing(Accumulator<RegOrder, Rec> &acc, int i)
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

        if constexpr (RegOrder >= 3)
        {
            acc.rho.uxuxux += moment_factor<MomentId::mxxx, MomentId::rho>(i);
            acc.rho.uxuxuy += moment_factor<MomentId::mxxy, MomentId::rho>(i);
            acc.rho.uxuyuy += moment_factor<MomentId::mxyy, MomentId::rho>(i);
            acc.rho.uyuyuy += moment_factor<MomentId::myyy, MomentId::rho>(i);
        }
    }

    template <>
    __device__ __forceinline__ void evaluate_outgoing(Accumulator<3, true> &acc, int i)
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

        acc.rho.uxuxux += (r::three * OMEGA - r::two) * moment_factor<MomentId::mxxx, MomentId::rho>(i);
        acc.rho.uxuxuy += (r::three * OMEGA - r::two) * moment_factor<MomentId::mxxy, MomentId::rho>(i);
        acc.rho.uxuyuy += (r::three * OMEGA - r::two) * moment_factor<MomentId::mxyy, MomentId::rho>(i);
        acc.rho.uyuyuy += (r::three * OMEGA - r::two) * moment_factor<MomentId::myyy, MomentId::rho>(i);

        acc.rho.uxmxx += (r::one - OMEGA) * r::three * moment_factor<MomentId::mxxx, MomentId::rho>(i);
        acc.rho.uymxx += (r::one - OMEGA) * moment_factor<MomentId::mxxy, MomentId::rho>(i);
        acc.rho.uxmxy += (r::one - OMEGA) * r::two * moment_factor<MomentId::mxxy, MomentId::rho>(i);
        acc.rho.uymxy += (r::one - OMEGA) * r::two * moment_factor<MomentId::mxyy, MomentId::rho>(i);
        acc.rho.uxmyy += (r::one - OMEGA) * moment_factor<MomentId::mxyy, MomentId::rho>(i);
        acc.rho.uymyy += (r::one - OMEGA) * r::three * moment_factor<MomentId::myyy, MomentId::rho>(i);
    }
}