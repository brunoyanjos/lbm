#pragma once

#include "core/types.cuh"

#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/factor.cuh"

#include "lbm/boundary/dirichlet/dirichlet_accumulator.cuh"

__device__ __forceinline__ void dirichlet_incoming_evaluation(DirichletAccumulator &acc,
                                                              const NodeMoments &M,
                                                              const real_t *__restrict__ pop, int i)
{
    acc.in.rho += pop[i];

    acc.in.mxx += pop[i] * hermite<MomentId::mxx>(i);
    acc.in.mxy += pop[i] * hermite<MomentId::mxy>(i);
    acc.in.myy += pop[i] * hermite<MomentId::myy>(i);

    acc.mxx.constant += moment_factor<MomentId::rho, MomentId::mxx>(i);
    acc.mxx.constant += M.ux * moment_factor<MomentId::ux, MomentId::mxx>(i);
    acc.mxx.constant += M.uy * moment_factor<MomentId::uy, MomentId::mxx>(i);

    acc.mxx.mxx += moment_factor<MomentId::mxx, MomentId::mxx>(i);
    acc.mxx.mxy += moment_factor<MomentId::mxy, MomentId::mxx>(i);
    acc.mxx.myy += moment_factor<MomentId::myy, MomentId::mxx>(i);

    acc.mxy.constant += moment_factor<MomentId::rho, MomentId::mxy>(i);
    acc.mxy.constant += M.ux * moment_factor<MomentId::ux, MomentId::mxy>(i);
    acc.mxy.constant += M.uy * moment_factor<MomentId::uy, MomentId::mxy>(i);

    acc.mxy.mxx += moment_factor<MomentId::mxx, MomentId::mxy>(i);
    acc.mxy.mxy += moment_factor<MomentId::mxy, MomentId::mxy>(i);
    acc.mxy.myy += moment_factor<MomentId::myy, MomentId::mxy>(i);

    acc.myy.constant += moment_factor<MomentId::rho, MomentId::myy>(i);
    acc.myy.constant += M.ux * moment_factor<MomentId::ux, MomentId::myy>(i);
    acc.myy.constant += M.uy * moment_factor<MomentId::uy, MomentId::myy>(i);

    acc.myy.mxx += moment_factor<MomentId::mxx, MomentId::myy>(i);
    acc.myy.mxy += moment_factor<MomentId::mxy, MomentId::myy>(i);
    acc.myy.myy += moment_factor<MomentId::myy, MomentId::myy>(i);
}