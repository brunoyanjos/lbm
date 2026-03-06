#pragma once
#include <cstddef>
#include "active_geometry.cuh"

constexpr real_t T_FINAL_STAR = static_cast<real_t>(30);
constexpr real_t SAVE_DT_STAR = static_cast<real_t>(1);

constexpr int N_STEPS = static_cast<int>(T_FINAL_STAR * (L_CHAR / U_MAX));
constexpr int SAVE_INTERVAL = static_cast<int>(SAVE_DT_STAR * (L_CHAR / U_MAX));
