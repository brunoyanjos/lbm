#pragma once

#include "../../core/types.cuh"
#include "lbm_state.cuh"
#include "../../core/physics.h"

template <int RegOrder, bool HighOrder>
__device__ __forceinline__ void store_next_state(const LBMStateFor<RegOrder, HighOrder> &S,
                                                 int n,
                                                 size_t idx,
                                                 const NodeMomentsFor<RegOrder, HighOrder> &M)
{
    S.d_rhoA[n][idx] = M.rhoA;
    S.d_rhoB[n][idx] = M.rhoB;

    S.d_ux[n][idx] = M.ux;
    S.d_uy[n][idx] = M.uy;
    S.d_mxx[n][idx] = M.mxx;
    S.d_mxy[n][idx] = M.mxy;
    S.d_myy[n][idx] = M.myy;
}
