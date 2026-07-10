#pragma once
#include "../../core/types.cuh"
#include "../../core/physics.h"

__device__ void moment_space_collision(NodeMoments &M)
{
    const real_t one_minus_omega = r::one - OMEGA;
    const real_t half_omega = r::half * OMEGA;

    M.mxx = one_minus_omega * M.mxx + half_omega * M.ux * M.ux;
    M.mxy = one_minus_omega * M.mxy + OMEGA * M.ux * M.uy;
    M.myy = one_minus_omega * M.myy + half_omega * M.uy * M.uy;
}