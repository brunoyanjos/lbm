#pragma once
#include <cstddef>
#include "geometry.h"
#include "physics.h"

constexpr int N_STEPS = real_t(250) * L_CHAR / U_MAX;
constexpr int SAVE_INTERVAL = L_CHAR / U_MAX;
