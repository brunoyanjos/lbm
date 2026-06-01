#pragma once

#include "core/physics.h"
#include "core/types.cuh"

#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/factor.cuh"

#include "lbm/boundary/dirichlet/accumulator.cuh"

namespace boundary::dirichlet
{
    template <int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ void evaluate_outgoing(Accumulator &acc,
                                                      const NodeMomentsFor<RegOrder, Rec, HighOrder> &M, int i)
    {
        acc.rho.constant += moment_factor<MomentId::rho, MomentId::rho>(i);
        acc.rho.constant += M.ux * moment_factor<MomentId::ux, MomentId::rho>(i);
        acc.rho.constant += M.uy * moment_factor<MomentId::uy, MomentId::rho>(i);
        acc.rho.constant += M.ux * M.ux * moment_factor<MomentId::mxx, MomentId::rho>(i);
        acc.rho.constant += OMEGA * M.ux * M.uy * moment_factor<MomentId::mxy, MomentId::rho>(i);
        acc.rho.constant += M.uy * M.uy * moment_factor<MomentId::myy, MomentId::rho>(i);

        if constexpr (RegOrder >= 3 && HighOrder)
        {
            acc.rho.constant += M.ux * M.ux * M.ux * moment_factor<MomentId::mxxx, MomentId::rho>(i);
            acc.rho.constant += M.uy * M.uy * M.uy * moment_factor<MomentId::myyy, MomentId::rho>(i);
        }

        if constexpr (RegOrder >= 3 && !Rec)
        {
            acc.rho.constant += M.ux * M.ux * M.uy * moment_factor<MomentId::mxxy, MomentId::rho>(i);
            acc.rho.constant += M.ux * M.uy * M.uy * moment_factor<MomentId::mxyy, MomentId::rho>(i);
        }

        if constexpr (RegOrder >= 3 && Rec)
        {
            acc.rho.constant += M.ux * M.ux * M.uy * (r::two * OMEGA - r::one) * moment_factor<MomentId::mxxy, MomentId::rho>(i);
            acc.rho.constant += M.ux * M.uy * M.uy * (r::two * OMEGA - r::one) * moment_factor<MomentId::mxyy, MomentId::rho>(i);
        }

        acc.rho.mxy += (r::one - OMEGA) * moment_factor<MomentId::mxy, MomentId::rho>(i);

        if constexpr (RegOrder >= 3 && Rec)
        {
            acc.rho.mxy += (r::one - OMEGA) * r::two * M.ux * moment_factor<MomentId::mxxy, MomentId::rho>(i);
            acc.rho.mxy += (r::one - OMEGA) * r::two * M.uy * moment_factor<MomentId::mxyy, MomentId::rho>(i);
        }
    }
}