#pragma once

#include "../../../core/types.cuh"
#include "moment_ids.cuh"
#include "factors.cuh"
#include "../../stencil_active.cuh"

__device__ __forceinline__ real_t linear_factor_value(MomentId id, const Factors &F)
{
    switch (id)
    {
    case MomentId::rho:
        return F.w;
    case MomentId::ux:
        return F.Hx;
    case MomentId::uy:
        return F.Hy;
    case MomentId::mxx:
        return F.Hxx;
    case MomentId::mxy:
        return F.Hxy;
    case MomentId::myy:
        return F.Hyy;

#if LBM_HAS_REG3_AXIAL
    case MomentId::mxxx:
        return F.Hxxx;
    case MomentId::myyy:
        return F.Hyyy;
#endif

#if LBM_HAS_REG3_CROSS
    case MomentId::mxxy:
        return F.Hxxy;
    case MomentId::mxyy:
        return F.Hxyy;
#endif

    default:
        return r::zero;
    }
}