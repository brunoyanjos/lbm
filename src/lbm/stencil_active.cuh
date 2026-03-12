#pragma once

#include "../core/lbm_stencil_select.cuh"

#if LBM_ACTIVE_STENCIL_D2Q9

#include "stencils/d2q9.cuh"
namespace Stencil = D2Q9;

#elif LBM_ACTIVE_STENCIL_D2V17

#include "stencils/d2v17.cuh"
namespace Stencil = D2V17;

#else
#error "No active stencil selected"
#endif