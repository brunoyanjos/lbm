#pragma once

#include "../../core/types.cuh"
#include "../stencil_active.cuh"
#include "scale_factor.cuh"
#include "moment_id.cuh"

__device__ __forceinline__ void scale_to_stored_basis(NodeMoments &M)
{
    M.ux *= scale_factor<MomentId::ux>();
    M.uy *= scale_factor<MomentId::uy>();
    M.mxx *= scale_factor<MomentId::mxx>();
    M.mxy *= scale_factor<MomentId::mxy>();
    M.myy *= scale_factor<MomentId::myy>();
}
