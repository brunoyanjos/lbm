#pragma once

// Escolha o stencil via macro de compilação:
//   -DLBM_STENCIL_D2Q9
//   -DLBM_STENCIL_D2V17
//   -DLBM_STENCIL_D2V37

#if defined(LBM_STENCIL_D2Q9)

#include "stencils/d2q9.cuh"
namespace Stencil = D2Q9;

#elif defined(LBM_STENCIL_D2V17)

#include "stencils/d2v17.cuh"
namespace Stencil = D2V17;

#elif defined(LBM_STENCIL_D2V37)

#include "stencils/d2v37.cuh"
namespace Stencil = D2V37;

#else

// default
#include "stencils/d2q9.cuh"
namespace Stencil = D2Q9;

#endif