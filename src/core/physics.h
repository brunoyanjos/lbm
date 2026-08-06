#pragma once
#include "types.cuh"
#include "geometry.h"
#include "math_utils.cuh"
#include "../lbm/stencil_active.cuh"

#ifndef LBM_RE
#define LBM_RE 10.0
#endif

constexpr real_t RE = static_cast<real_t>(LBM_RE);

constexpr real_t GRAVITY_Y = -static_cast<real_t>(0.000008);
constexpr real_t GRAVITY_X = r::zero;
constexpr real_t L_CHAR = static_cast<real_t>(NX - 1);

constexpr real_t U = constexpr_sqrt(constexpr_sqrt(GRAVITY_Y * GRAVITY_Y + GRAVITY_X * GRAVITY_X) * L_CHAR);

constexpr real_t VISC = U * L_CHAR / RE;
constexpr real_t TAU = static_cast<real_t>(0.5) + Stencil::as2 * VISC;
constexpr real_t OMEGA = static_cast<real_t>(1.0) / TAU;

constexpr real_t RHOA_0 = static_cast<real_t>(1.0);
constexpr real_t RHOB_0 = static_cast<real_t>(1.0);

constexpr real_t BETA = static_cast<real_t>(0.7);
