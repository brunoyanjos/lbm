#pragma once

#include "../../core/local_domain.cuh"
#include "../../core/simulation_config.h"
#include "../../core/types.cuh"
#include "lbm/stencil_active.cuh"

struct SecondOrderLocalState
{
    LocalDomain domain;

    // host (local buffer, including halo rows)
    real_t *h_rho, *h_ux, *h_uy, *h_mxx, *h_mxy, *h_myy;

    // device ping-pong (local buffer, including halo rows)
    real_t *d_rho[2], *d_ux[2], *d_uy[2];
    real_t *d_mxx[2], *d_mxy[2], *d_myy[2];

    int cur;
    size_t N;
    size_t bytes_field;
};

template <int Order, bool Rec, bool HighOrder>
struct LBMLocalStateFor : SecondOrderLocalState
{
};

template <>
struct LBMLocalStateFor<3, false, false> : SecondOrderLocalState
{
    real_t *h_mxxy, *h_mxyy;
    real_t *d_mxxy[2], *d_mxyy[2];
};

template <>
struct LBMLocalStateFor<3, false, true> : SecondOrderLocalState
{
    real_t *h_mxxx, *h_mxxy, *h_mxyy, *h_myyy;
    real_t *d_mxxx[2], *d_mxxy[2], *d_mxyy[2], *d_myyy[2];
};

using LBMLocalState = LBMLocalStateFor<REG_ORDER, USE_RECURRENCE, Stencil::high_order>;

[[nodiscard]] __host__ LBMLocalState lbm_allocate_local_state(const LocalDomain &domain);

__host__ void lbm_free_local_state(LBMLocalState &S);
