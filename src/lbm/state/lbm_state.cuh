#pragma once

#include "../../core/types.cuh"
#include "../../core/simulation_config.h"
#include "../../app/cuda_config.cuh"
#include "lbm/stencil_active.cuh"

struct SecondOrderState
{
    // host (single buffer for output)
    real_t *h_rhoA, *h_uxA, *h_uyA, *h_mxxA, *h_mxyA, *h_myyA;
    real_t *h_rhoB, *h_uxB, *h_uyB, *h_mxxB, *h_mxyB, *h_myyB;

    // device ping-pong
    real_t *d_rhoA[2], *d_uxA[2], *d_uyA[2];
    real_t *d_mxxA[2], *d_mxyA[2], *d_myyA[2];

    real_t *d_rhoB[2], *d_uxB[2], *d_uyB[2];
    real_t *d_mxxB[2], *d_mxyB[2], *d_myyB[2];

    int cur;
    size_t N;
    size_t bytes_field;
};

template <int Order, bool Rec, bool HighOrder>
struct LBMStateFor : SecondOrderState
{
};

template <>
struct LBMStateFor<3, false, false> : SecondOrderState
{
    real_t *h_mxxy, *h_mxyy;
    real_t *d_mxxy[2], *d_mxyy[2];
};

template <>
struct LBMStateFor<3, false, true> : SecondOrderState
{
    real_t *h_mxxx, *h_mxxy, *h_mxyy, *h_myyy;
    real_t *d_mxxx[2], *d_mxxy[2], *d_mxyy[2], *d_myyy[2];
};

using LBMState = LBMStateFor<REG_ORDER, USE_RECURRENCE, Stencil::high_order>;

[[nodiscard]] __host__ LBMState lbm_allocate_state();

__host__ void lbm_free_state(LBMState &S);
