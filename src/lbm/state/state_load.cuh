#pragma once

#include "../../core/physics.h"
#include "lbm/moment/moment_id.cuh"
#include "lbm/moment/scale_factor.cuh"
#include "lbm/state/lbm_state.cuh"

template <int RegOrder, bool Rec, bool HighOrder>
__device__ __forceinline__ void load_state_moments(const LBMStateFor<RegOrder, Rec, HighOrder> &S,
                                                   int c,
                                                   size_t idx,
                                                   NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
{
    M.rho = S.d_rho[c][idx] + RHOA_0;
    M.ux = S.d_ux[c][idx] * inv_scale_factor<MomentId::ux>();
    M.uy = S.d_uy[c][idx] * inv_scale_factor<MomentId::uy>();
    M.mxx = S.d_mxx[c][idx] * inv_scale_factor<MomentId::mxx>();
    M.mxy = S.d_mxy[c][idx] * inv_scale_factor<MomentId::mxy>();
    M.myy = S.d_myy[c][idx] * inv_scale_factor<MomentId::myy>();
}

template <bool HighOrder>
__device__ __forceinline__ void load_state_moments(const LBMStateFor<3, false, HighOrder> &S,
                                                   int c,
                                                   size_t idx,
                                                   NodeMomentsFor<3, false, HighOrder> &M)
{
    M.rho = S.d_rho[c][idx] + RHOA_0;
    M.ux = S.d_ux[c][idx] * inv_scale_factor<MomentId::ux>();
    M.uy = S.d_uy[c][idx] * inv_scale_factor<MomentId::uy>();
    M.mxx = S.d_mxx[c][idx] * inv_scale_factor<MomentId::mxx>();
    M.mxy = S.d_mxy[c][idx] * inv_scale_factor<MomentId::mxy>();
    M.myy = S.d_myy[c][idx] * inv_scale_factor<MomentId::myy>();
    M.mxxy = S.d_mxxy[c][idx] * inv_scale_factor<MomentId::mxxy>();
    M.mxyy = S.d_mxyy[c][idx] * inv_scale_factor<MomentId::mxyy>();

    if constexpr (HighOrder)
    {
        M.mxxx = S.d_mxxx[c][idx] * inv_scale_factor<MomentId::mxxx>();
        M.myyy = S.d_myyy[c][idx] * inv_scale_factor<MomentId::myyy>();
    }
}
