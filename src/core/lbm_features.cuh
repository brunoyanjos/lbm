#pragma once

#ifndef LBM_REG_ORDER
#define LBM_REG_ORDER 3
#endif

static_assert(LBM_REG_ORDER == 2 || LBM_REG_ORDER == 3,
              "LBM_REG_ORDER must be 2 or 3");

#if defined(LBM_STENCIL_D2Q9)
#define LBM_BUILD_STENCIL_D2Q9 1
#define LBM_BUILD_STENCIL_D2V17 0
#elif defined(LBM_STENCIL_D2V17)
#define LBM_BUILD_STENCIL_D2Q9 0
#define LBM_BUILD_STENCIL_D2V17 1
#else
// fallback para editor / default local
#define LBM_BUILD_STENCIL_D2Q9 1
#define LBM_BUILD_STENCIL_D2V17 0
#endif

#define LBM_HAS_MXXY_MXYY 1
#define LBM_HAS_MXXX_MYYY LBM_BUILD_STENCIL_D2V17

#if LBM_REG_ORDER >= 3
#define LBM_HAS_REG3_CROSS LBM_HAS_MXXY_MXYY
#define LBM_HAS_REG3_AXIAL LBM_HAS_MXXX_MYYY
#else
#define LBM_HAS_REG3_CROSS 0
#define LBM_HAS_REG3_AXIAL 0
#endif