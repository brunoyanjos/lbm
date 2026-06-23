#pragma once

#include "../../core/types.cuh"
#include "lbm_state.cuh"
#include "../../core/physics.h"

template <int RegOrder, bool Rec, bool HighOrder>
__device__ __forceinline__ void store_next_state(const LBMStateFor<RegOrder, Rec, HighOrder> &S,
                                                 int n,
                                                 size_t idx,
                                                 const NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
{
    S.d_rho[n][idx] = M.rho;// - RHO_0;
    S.d_ux[n][idx] = M.ux;
    S.d_uy[n][idx] = M.uy;
    S.d_mxx[n][idx] = M.mxx;
    S.d_mxy[n][idx] = M.mxy;
    S.d_myy[n][idx] = M.myy;
}

template <bool HighOrder>
__device__ __forceinline__ void store_next_state(const LBMStateFor<3, false, HighOrder> &S,
                                                 int n,
                                                 size_t idx,
                                                 const NodeMomentsFor<3, false, HighOrder> &M)
{
    S.d_rho[n][idx] = M.rho;// - RHO_0;
    S.d_ux[n][idx] = M.ux;
    S.d_uy[n][idx] = M.uy;
    S.d_mxx[n][idx] = M.mxx;
    S.d_mxy[n][idx] = M.mxy;
    S.d_myy[n][idx] = M.myy;
    S.d_mxxy[n][idx] = M.mxxy;
    S.d_mxyy[n][idx] = M.mxyy;

    if constexpr (HighOrder)
    {
        S.d_mxxx[n][idx] = M.mxxx;
        S.d_myyy[n][idx] = M.myyy;
    }
}
