#pragma once

#include "core/types.cuh"

#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/accumulator.cuh"
#include "lbm/boundary/common/factor.cuh"

#include "lbm/boundary/dirichlet/dirichlet_accumulator.cuh"

__device__ __forceinline__ void dirichlet_incoming_evaluation(DirichletAccumulator &acc,
                                                              const NodeMoments &M,
                                                              const real_t *__restrict__ pop, int i)
{
    acc.in.rho += pop[i];
    acc.in.mxy += pop[i] * hermite<MomentId::mxy>(i);

    acc.mxy.constant += moment_factor<MomentId::rho, MomentId::mxy>(i);
    acc.mxy.constant += M.ux * moment_factor<MomentId::ux, MomentId::mxy>(i);
    acc.mxy.constant += M.uy * moment_factor<MomentId::uy, MomentId::mxy>(i);
    acc.mxy.constant += M.ux * M.ux * moment_factor<MomentId::mxx, MomentId::mxy>(i);
    acc.mxy.constant += M.uy * M.uy * moment_factor<MomentId::myy, MomentId::mxy>(i);

    acc.mxy.mxy += moment_factor<MomentId::mxy, MomentId::mxy>(i);
}