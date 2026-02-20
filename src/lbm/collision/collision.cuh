#pragma once
#include "../../core/types.cuh"
#include "../../core/physics.h"

__device__ void moment_space_collision(const real_t &ux, const real_t &uy,
                                       real_t &mxx, real_t &mxy, real_t &myy)
{
    const real_t one_minus_omega = static_cast<real_t>(1) - OMEGA;
    const real_t half_omega = static_cast<real_t>(0.5) * OMEGA;

    mxx = one_minus_omega * mxx + half_omega * ux * ux;
    mxy = one_minus_omega * mxy + OMEGA * ux * uy;
    myy = one_minus_omega * myy + half_omega * uy * uy;
}