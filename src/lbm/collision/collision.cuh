#pragma once
#include "../../core/types.cuh"
#include "../../core/physics.h"

template <int RegOrder, bool Rec, bool HighOrder>
__device__ void moment_space_collision(NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
{
    const real_t one_minus_omega = r::one - OMEGA;
    const real_t half_omega = r::half * OMEGA;

    M.mxx = one_minus_omega * M.mxx + half_omega * M.ux * M.ux;
    M.mxy = one_minus_omega * M.mxy + OMEGA * M.ux * M.uy;
    M.myy = one_minus_omega * M.myy + half_omega * M.uy * M.uy;
}

template <bool HighOrder>
__device__ void moment_space_collision(NodeMomentsFor<3, false, HighOrder> &M)
{
    const real_t one_minus_omega = r::one - OMEGA;
    const real_t half_omega = r::half * OMEGA;

    M.mxx = one_minus_omega * M.mxx + half_omega * M.ux * M.ux;
    M.mxy = one_minus_omega * M.mxy + OMEGA * M.ux * M.uy;
    M.myy = one_minus_omega * M.myy + half_omega * M.uy * M.uy;

    M.mxxy = one_minus_omega * M.mxxy + half_omega * M.ux * M.ux * M.uy;
    M.mxyy = one_minus_omega * M.mxyy + half_omega * M.ux * M.uy * M.uy;

    if constexpr (HighOrder)
    {
        const real_t sixth_omega = r::sixth * OMEGA;
        M.mxxx = one_minus_omega * M.mxxx + sixth_omega * M.ux * M.ux * M.ux;
        M.myyy = one_minus_omega * M.myyy + sixth_omega * M.uy * M.uy * M.uy;
    }
}
