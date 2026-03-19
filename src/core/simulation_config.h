#pragma once
#include <cstddef>
#include "geometry.h"
#include "physics.h"

constexpr int N_STEPS = real_t(1000) * NX / U_LID;
constexpr int SAVE_INTERVAL = NX / U_LID;
