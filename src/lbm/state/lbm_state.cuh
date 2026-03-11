#pragma once

#include "../../core/types.cuh"
#include "../../app/cuda_config.cuh"

struct LBMState
{
    // host (single buffer for output)
    real_t *h_rho, *h_ux, *h_uy, *h_mxx, *h_mxy, *h_myy;

    real_t *global_rho, *global_ux, *global_uy;
    real_t *global_mxx, *global_mxy, *global_myy;

    // device ping-pong
    real_t *d_rho, *d_ux, *d_uy;
    real_t *d_mxx, *d_mxy, *d_myy;

    size_t N;
    size_t bytes_field;
};

[[nodiscard]] __host__ LBMState lbm_allocate_state();

__host__ void lbm_free_state(LBMState &S);