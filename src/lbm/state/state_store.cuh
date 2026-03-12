#pragma once

#include "../../core/types.cuh"
#include "lbm_state.cuh"
#include "../../geometries/active_geometry.cuh"

__device__ __forceinline__ void store_next_state(const LBMState &S,
                                                 int n,
                                                 size_t idx,
                                                 real_t rho,
                                                 real_t ux, real_t uy,
                                                 real_t mxx, real_t mxy, real_t myy)
{
    S.d_rho[n][idx] = rho - Geometry::RHO_0;
    S.d_ux[n][idx] = ux;
    S.d_uy[n][idx] = uy;
    S.d2.mxx[n][idx] = mxx;
    S.d2.mxy[n][idx] = mxy;
    S.d2.myy[n][idx] = myy;
}
