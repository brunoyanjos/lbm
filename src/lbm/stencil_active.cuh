#pragma once

// Escolha o stencil via macro de compilação:
//   -DLBM_STENCIL_D2Q9  ou  -DLBM_STENCIL_D2V17
#if defined(LBM_STENCIL_D2Q9)
#include "stencils/d2q9.cuh"
namespace Stencil = D2Q9;
#elif defined(LBM_STENCIL_D2V17)
#include "stencils/d2v17.cuh"
namespace Stencil = D2V17;
#else
// default
#include "stencils/d2q9.cuh"
namespace Stencil = D2Q9;
#endif
