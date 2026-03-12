#pragma once

#include "../../../core/types.cuh"
#include "nonlinear_ids.cuh"
#include "factors.cuh"
#include "../../stencil_active.cuh"

__device__ __forceinline__ real_t nonlinear_factor_value(NonlinearId id, const Factors &F)
{
    switch (id)
    {
    case NonlinearId::uxux:
        return F.Hxx;
    case NonlinearId::uxuy:
        return F.Hxy;
    case NonlinearId::uyuy:
        return F.Hyy;

#if LBM_HAS_REG3_AXIAL
    case NonlinearId::uxuxux:
        return F.Hxxx;
    case NonlinearId::uyuyuy:
        return F.Hyyy;
#endif

#if LBM_HAS_REG3_CROSS
    case NonlinearId::uxuxuy:
        return F.Hxxy;
    case NonlinearId::uxuyuy:
        return F.Hxyy;
#endif

    default:
        return r::zero;
    }
}