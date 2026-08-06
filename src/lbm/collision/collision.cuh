#pragma once
#include "lbm/moment/scale_factor.cuh"
#include "../../core/types.cuh"
#include "../../core/physics.h"

template <int RegOrder, bool Rec, bool HighOrder>
__device__ void moment_space_collision(NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
{
    const real_t ux_col = M.ux - r::half * GRAVITY_X;
    const real_t uy_col = M.uy - r::half * GRAVITY_Y;

    const real_t mxx_col = (r::one - OMEGA) * M.mxx + OMEGA * M.ux * M.ux - GRAVITY_X * M.ux;
    const real_t mxy_col = (r::one - OMEGA) * M.mxy + OMEGA * M.ux * M.uy - r::half * (GRAVITY_Y * M.ux + GRAVITY_X * M.uy);
    const real_t myy_col = (r::one - OMEGA) * M.myy + OMEGA * M.uy * M.uy - GRAVITY_Y * M.uy;

    M.ux = ux_col;
    M.uy = uy_col;

    M.mxx = mxx_col;
    M.mxy = mxy_col;
    M.myy = myy_col;
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
