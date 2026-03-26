#pragma once

#include "core/physics.h"
#include "core/types.cuh"

#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/factor.cuh"

#include "lbm/boundary/dirichlet/dirichlet_accumulator.cuh"

__device__ __forceinline__ void dirichlet_outgoing_evaluation(DirichletAccumulator &acc,
                                                              const NodeMoments &M, int i)
{
    acc.rho.constant += moment_factor<MomentId::rho, MomentId::rho>(i);
    acc.rho.constant += M.ux * moment_factor<MomentId::ux, MomentId::rho>(i);
    acc.rho.constant += M.uy * moment_factor<MomentId::uy, MomentId::rho>(i);
    acc.rho.constant += OMEGA * M.ux * M.ux * moment_factor<MomentId::mxx, MomentId::rho>(i);
    acc.rho.constant += OMEGA * M.ux * M.uy * moment_factor<MomentId::mxy, MomentId::rho>(i);
    acc.rho.constant += OMEGA * M.uy * M.uy * moment_factor<MomentId::myy, MomentId::rho>(i);

    acc.rho.mxx += (r::one - OMEGA) * moment_factor<MomentId::mxx, MomentId::rho>(i);
    acc.rho.mxy += (r::one - OMEGA) * moment_factor<MomentId::mxy, MomentId::rho>(i);
    acc.rho.myy += (r::one - OMEGA) * moment_factor<MomentId::myy, MomentId::rho>(i);
}