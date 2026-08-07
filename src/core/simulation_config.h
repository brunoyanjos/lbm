#pragma once
#include <cstddef>
#include "geometry.h"
#include "physics.h"

#ifndef LBM_T_STAR_END
#define LBM_T_STAR_END 1000
#endif

#ifndef LBM_N_STEPS
#define LBM_N_STEPS (real_t(LBM_T_STAR_END) * NX / U)
#endif

#ifndef LBM_SAVE_INTERVAL
#define LBM_SAVE_INTERVAL (NX / U)
#endif

#ifndef LBM_VTI_SAVE_INTERVAL
#define LBM_VTI_SAVE_INTERVAL (100 * LBM_SAVE_INTERVAL)
#endif

#ifndef LBM_AVG_START_T_STAR
#define LBM_AVG_START_T_STAR 0
#endif

#ifndef LBM_REG_ORDER
#define LBM_REG_ORDER 2
#endif

#if LBM_REG_ORDER != 2 && LBM_REG_ORDER != 3
#error "LBM_REG_ORDER must be 2 or 3"
#endif

#ifndef LBM_USE_RECURRENCE
#define LBM_USE_RECURRENCE 0
#endif

#if LBM_USE_RECURRENCE != 0 && LBM_USE_RECURRENCE != 1
#error "LBM_USE_RECURRENCE must be 0 or 1"
#endif

#if LBM_USE_RECURRENCE == 1 && LBM_REG_ORDER != 3
#error "LBM_USE_RECURRENCE requires LBM_REG_ORDER == 3"
#endif

#ifndef LBM_USE_SYMBOLIC_BOUNDARY
#define LBM_USE_SYMBOLIC_BOUNDARY 0
#endif

#if LBM_USE_SYMBOLIC_BOUNDARY != 0 && LBM_USE_SYMBOLIC_BOUNDARY != 1
#error "LBM_USE_SYMBOLIC_BOUNDARY must be 0 or 1"
#endif

constexpr int N_STEPS = 200000;
constexpr int SAVE_INTERVAL = 50000;
constexpr int VTI_SAVE_INTERVAL = 50000;
constexpr int AVG_START_STEP = 0;
constexpr int REG_ORDER = LBM_REG_ORDER;
constexpr bool USE_RECURRENCE = LBM_USE_RECURRENCE != 0;
constexpr bool USE_SYMBOLIC_BOUNDARY = LBM_USE_SYMBOLIC_BOUNDARY != 0;
