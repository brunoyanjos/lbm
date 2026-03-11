#pragma once

#if defined(LBM_STENCIL_D2Q9)
#include "stencils/d2q9.cuh"
namespace Stencil = D2Q9;
#define LBM_ACTIVE_STENCIL_D2Q9 1
#define LBM_ACTIVE_STENCIL_D2V17 0
#elif defined(LBM_STENCIL_D2V17)
#include "stencils/d2v17.cuh"
namespace Stencil = D2V17;
#define LBM_ACTIVE_STENCIL_D2Q9 0
#define LBM_ACTIVE_STENCIL_D2V17 1
#else
#include "stencils/d2q9.cuh"
namespace Stencil = D2Q9;
#define LBM_ACTIVE_STENCIL_D2Q9 1
#define LBM_ACTIVE_STENCIL_D2V17 0
#endif