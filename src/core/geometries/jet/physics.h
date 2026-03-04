#pragma once
#include "../../types.cuh"
#include "geometry.h"
#include "../../../lbm/stencil_active.cuh"

#ifndef LBM_RE
#define LBM_RE 100.0
#endif

constexpr real_t RE = static_cast<real_t>(LBM_RE);

constexpr real_t U_MAX = static_cast<real_t>(0.0256);
constexpr real_t L_CHAR = static_cast<real_t>(JET_H);

constexpr real_t VISC = (U_MAX * L_CHAR) / RE;

constexpr real_t TAU = static_cast<real_t>(0.5) + Stencil::as2 * VISC;
constexpr real_t OMEGA = static_cast<real_t>(1.0) / TAU;

constexpr real_t RHO_0 = static_cast<real_t>(1.0);