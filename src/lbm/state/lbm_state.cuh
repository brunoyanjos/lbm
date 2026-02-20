#pragma once

#include "../../core/types.cuh"
#include "../../app/cuda_config.cuh"

struct LBMState
{
    // host (single buffer for output)
    real_t *h_rho, *h_ux, *h_uy, *h_mxx, *h_mxy, *h_myy;

    // device ping-pong
    real_t *d_rho[2], *d_ux[2], *d_uy[2];
    real_t *d_mxx[2], *d_mxy[2], *d_myy[2];

    int cur; // 0 or 1
    size_t N;
    size_t bytes_field;
};

[[nodiscard]] __host__ LBMState lbm_allocate_state();

__host__ void lbm_free_state(LBMState &S);