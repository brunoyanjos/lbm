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
    S.d_rho[idx] = rho - Geometry::RHO_0;
    S.d_ux[idx] = ux;
    S.d_uy[idx] = uy;
    S.d_mxx[idx] = mxx;
    S.d_mxy[idx] = mxy;
    S.d_myy[idx] = myy;
}
