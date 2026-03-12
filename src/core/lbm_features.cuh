#pragma once

#include "lbm_stencil_select.cuh"

#ifndef LBM_REG_ORDER
#define LBM_REG_ORDER 2
#endif

static_assert(LBM_REG_ORDER == 2 || LBM_REG_ORDER == 3,
              "LBM_REG_ORDER must be 2 or 3");

#define LBM_HAS_MXXY_MXYY 1
#define LBM_HAS_MXXX_MYYY LBM_ACTIVE_STENCIL_D2V17

#if LBM_REG_ORDER >= 3
#define LBM_HAS_REG3_CROSS LBM_HAS_MXXY_MXYY
#define LBM_HAS_REG3_AXIAL LBM_HAS_MXXX_MYYY
#else
#define LBM_HAS_REG3_CROSS 0
#define LBM_HAS_REG3_AXIAL 0
#endif