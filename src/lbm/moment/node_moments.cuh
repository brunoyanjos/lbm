#pragma once

#include "core/types.cuh"
#include "core/simulation_config.h"
#include "lbm/stencil_active.cuh"

struct SecondOrderMoments
{
    real_t rhoA, rhoB;
    real_t ux;
    real_t uy;
    real_t mxx;
    real_t mxy;
    real_t myy;
};

template <int Order, bool HighOrder>
struct NodeMomentsFor : SecondOrderMoments
{
};

template <>
struct NodeMomentsFor<3, false> : SecondOrderMoments
{
    real_t mxxy;
    real_t mxyy;
};

template <>
struct NodeMomentsFor<3, true> : SecondOrderMoments
{
    real_t mxxx;
    real_t mxxy;
    real_t mxyy;
    real_t myyy;
};

using NodeMoments = NodeMomentsFor<REG_ORDER, Stencil::high_order>;
