#pragma once

#include "../../core/types.cuh"
#include "lbm_state.cuh"
#include "../../core/physics.h"

__device__ __forceinline__ void store_next_state(const LBMState &S,
                                                 int n,
                                                 size_t idx,
                                                 const NodeMoments &M)
{
    S.d_rho[n][idx] = M.rho - RHO_0;
    S.d_ux[n][idx] = M.ux;
    S.d_uy[n][idx] = M.uy;
    S.d_mxx[n][idx] = M.mxx;
    S.d_mxy[n][idx] = M.mxy;
    S.d_myy[n][idx] = M.myy;
}
