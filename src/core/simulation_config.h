#pragma once
#include <cstddef>
#include "geometry.h"
#include "physics.h"

constexpr real_t T_FINAL_STAR = static_cast<real_t>(10);
constexpr real_t SAVE_DT_STAR = static_cast<real_t>(1);

constexpr int N_STEPS = static_cast<int>(T_FINAL_STAR * (L_CHAR / U_WALL));
constexpr int SAVE_INTERVAL = static_cast<int>(SAVE_DT_STAR * (L_CHAR / U_WALL));
