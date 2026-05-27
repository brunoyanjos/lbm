#pragma once
#include <cstddef>
#include "geometry.h"
#include "physics.h"

#ifndef LBM_T_STAR_END
#define LBM_T_STAR_END 1000
#endif

#ifndef LBM_N_STEPS
#define LBM_N_STEPS (real_t(LBM_T_STAR_END) * NX / U_LID)
#endif

#ifndef LBM_SAVE_INTERVAL
#define LBM_SAVE_INTERVAL (NX / U_LID)
#endif

#ifndef LBM_VTI_SAVE_INTERVAL
#define LBM_VTI_SAVE_INTERVAL (100 * LBM_SAVE_INTERVAL)
#endif

#ifndef LBM_AVG_START_T_STAR
#define LBM_AVG_START_T_STAR 0
#endif

constexpr int N_STEPS = LBM_N_STEPS;
constexpr int SAVE_INTERVAL = LBM_SAVE_INTERVAL;
constexpr int VTI_SAVE_INTERVAL = LBM_VTI_SAVE_INTERVAL;
constexpr int AVG_START_STEP = real_t(LBM_AVG_START_T_STAR) * NX / U_LID;
