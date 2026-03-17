#pragma once

#include "core/types.cuh"
#include "core/geometry.h"
#include "core/physics.h"

#include "lbm/moment/node_moments.cuh"

__device__ __forceinline__ void bc_velocity(NodeMoments &M, const int &x, const int &y)
{
    M.ux = r::zero;
    M.uy = r::zero;

    if (y == NY - 1)
        M.ux = U_LID;
}