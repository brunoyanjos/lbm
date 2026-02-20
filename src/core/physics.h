#pragma once
#include "types.cuh"
#include "geometry.h"

constexpr real_t RE = static_cast<real_t>(1000.0);
constexpr real_t U_LID = static_cast<real_t>(0.0256);

constexpr real_t L_CHAR = static_cast<real_t>(NX - 1);

constexpr real_t VISC = U_LID * L_CHAR / RE;
constexpr real_t TAU = static_cast<real_t>(0.5) + static_cast<real_t>(3.0) * VISC;
constexpr real_t OMEGA = static_cast<real_t>(1.0) / TAU;

constexpr real_t RHO_0 = static_cast<real_t>(1.0);