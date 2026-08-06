#pragma once

#include "../../core/types.cuh"
#include "../../core/simulation_config.h"
#include "../stencil_active.cuh"
#include "scale_factor.cuh"
#include "moment_id.cuh"

template <int RegOrder, bool HighOrder>
__device__ __forceinline__ void scale_to_stored_basis(NodeMomentsFor<RegOrder, HighOrder> &M)
{
    M.uxA *= scale_factor<MomentId::ux>();
    M.uyA *= scale_factor<MomentId::uy>();
    M.mxxA *= scale_factor<MomentId::mxx>();
    M.mxyA *= scale_factor<MomentId::mxy>();
    M.myyA *= scale_factor<MomentId::myy>();
}