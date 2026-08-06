#pragma once

#include "../../core/types.cuh"
#include "../../core/simulation_config.h"
#include "../stencil_active.cuh"
#include "scale_factor.cuh"
#include "moment_id.cuh"

template <int RegOrder, bool Rec, bool HighOrder>
__device__ __forceinline__ void scale_to_stored_basis(NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
{
    M.uxA *= scale_factor<MomentId::ux>();
    M.uyA *= scale_factor<MomentId::uy>();
    M.mxxA *= scale_factor<MomentId::mxx>();
    M.mxyA *= scale_factor<MomentId::mxy>();
    M.myyA *= scale_factor<MomentId::myy>();
}

template <bool HighOrder>
__device__ __forceinline__ void scale_to_stored_basis(NodeMomentsFor<3, false, HighOrder> &M)
{
    M.ux *= scale_factor<MomentId::ux>();
    M.uy *= scale_factor<MomentId::uy>();
    M.mxx *= scale_factor<MomentId::mxx>();
    M.mxy *= scale_factor<MomentId::mxy>();
    M.myy *= scale_factor<MomentId::myy>();

    M.mxxy *= scale_factor<MomentId::mxxy>();
    M.mxyy *= scale_factor<MomentId::mxyy>();
    if constexpr (HighOrder)
    {
        M.mxxx *= scale_factor<MomentId::mxxx>();
        M.myyy *= scale_factor<MomentId::myyy>();
    }
}
