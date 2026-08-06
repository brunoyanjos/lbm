#pragma once

#include "../../core/physics.h"
#include "lbm/moment/moment_id.cuh"
#include "lbm/moment/scale_factor.cuh"
#include "lbm/state/lbm_state.cuh"

template <int RegOrder, bool HighOrder>
__device__ __forceinline__ void load_state_moments(const LBMStateFor<RegOrder, HighOrder> &S,
                                                   int c,
                                                   size_t idx,
                                                   NodeMomentsFor<RegOrder, HighOrder> &M)
{
    M.rho = S.d_rho[c][idx] + RHOA_0;
    M.ux = S.d_ux[c][idx] * inv_scale_factor<MomentId::ux>();
    M.uy = S.d_uy[c][idx] * inv_scale_factor<MomentId::uy>();
    M.mxx = S.d_mxx[c][idx] * inv_scale_factor<MomentId::mxx>();
    M.mxy = S.d_mxy[c][idx] * inv_scale_factor<MomentId::mxy>();
    M.myy = S.d_myy[c][idx] * inv_scale_factor<MomentId::myy>();
}
