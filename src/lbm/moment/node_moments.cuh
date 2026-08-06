#pragma once

#include "core/types.cuh"
#include "core/simulation_config.h"
#include "lbm/stencil_active.cuh"

struct SecondOrderMoments
{
    real_t rhoA, rhoB;
    real_t uxA, uxB;
    real_t uyA, uyB;
    real_t mxxA, mxxB;
    real_t mxyA, mxyB;
    real_t myyA, myyB;
};

template <int Order, bool Rec, bool HighOrder>
struct NodeMomentsFor : SecondOrderMoments
{
};

template <>
struct NodeMomentsFor<3, false, false> : SecondOrderMoments
{
    real_t mxxy;
    real_t mxyy;
};

template <>
struct NodeMomentsFor<3, false, true> : SecondOrderMoments
{
    real_t mxxx;
    real_t mxxy;
    real_t mxyy;
    real_t myyy;
};

using NodeMoments = NodeMomentsFor<REG_ORDER, USE_RECURRENCE, Stencil::high_order>;
