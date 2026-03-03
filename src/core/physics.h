#pragma once
#include "types.cuh"
#include "geometry.h"
#include "../lbm/stencil_active.cuh"

#ifndef LBM_RE
#define LBM_RE 7500.0
#endif

constexpr real_t RE = static_cast<real_t>(LBM_RE);
constexpr real_t U_WALL = static_cast<real_t>(0.0256);

constexpr real_t R_IN = static_cast<real_t>(NX) * static_cast<real_t>(0.25);
constexpr real_t R_OUT = static_cast<real_t>(NX - 1) * static_cast<real_t>(0.5);
constexpr real_t L_CHAR = (R_OUT - R_IN);

constexpr real_t VISC = (U_WALL * L_CHAR) / RE;

constexpr real_t TAU = static_cast<real_t>(0.5) + Stencil::as2 * VISC;
constexpr real_t OMEGA = static_cast<real_t>(1.0) / TAU;

constexpr real_t RHO_0 = static_cast<real_t>(1.0);