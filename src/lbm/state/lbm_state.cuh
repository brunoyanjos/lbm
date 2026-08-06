#pragma once

#include "../../core/types.cuh"
#include "../../core/simulation_config.h"
#include "../../app/cuda_config.cuh"
#include "lbm/interface/normal.cuh"
#include "lbm/stencil_active.cuh"

struct SecondOrderState
{
    // host (single buffer for output)
    real_t *h_rhoA, *h_rhoB, *h_ux, *h_uy, *h_mxx, *h_mxy, *h_myy;

    // device ping-pong
    real_t *d_rhoA[2], *d_rhoB[2], *d_ux[2], *d_uy[2];
    real_t *d_mxx[2], *d_mxy[2], *d_myy[2];

    Normal *d_n[2];

    int cur;
    size_t N;
    size_t bytes_field;
};

template <int Order, bool HighOrder>
struct LBMStateFor : SecondOrderState
{
};

template <>
struct LBMStateFor<3, false> : SecondOrderState
{
    real_t *h_mxxy, *h_mxyy;
    real_t *d_mxxy[2], *d_mxyy[2];
};

template <>
struct LBMStateFor<3, true> : SecondOrderState
{
    real_t *h_mxxx, *h_mxxy, *h_mxyy, *h_myyy;
    real_t *d_mxxx[2], *d_mxxy[2], *d_mxyy[2], *d_myyy[2];
};

using LBMState = LBMStateFor<REG_ORDER, Stencil::high_order>;

[[nodiscard]] __host__ LBMState lbm_allocate_state();

__host__ void lbm_free_state(LBMState &S);
