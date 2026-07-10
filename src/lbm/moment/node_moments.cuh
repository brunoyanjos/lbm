#pragma once

#include "core/types.cuh"
#include "core/simulation_config.h"
#include "lbm/stencil_active.cuh"

struct NodeMoments
{
    real_t rho;
    real_t ux;
    real_t uy;
    real_t mxx;
    real_t mxy;
    real_t myy;
};
