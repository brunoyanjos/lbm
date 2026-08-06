#pragma once

#include "core/types.cuh"

#include "lbm/stencil_active.cuh"
#include "lbm/moment/moment_id.cuh"
#include "lbm/meta/always_false.cuh"

template <auto Id>
__host__ __device__ __forceinline__
    real_t
    hermite(real_t cx, real_t cy)
{
    if constexpr (Id == MomentId::rho)
        return r::one;
    else if constexpr (Id == MomentId::ux)
        return cx;
    else if constexpr (Id == MomentId::uy)
        return cy;
    else if constexpr (Id == MomentId::mxx)
        return cx * cx - Stencil::cs2;
    else if constexpr (Id == MomentId::mxy)
        return cx * cy;
    else if constexpr (Id == MomentId::myy)
        return cy * cy - Stencil::cs2;
    else if constexpr (Id == MomentId::mxxx)
        return cx * cx * cx - r::three * Stencil::cs2 * cx;
    else if constexpr (Id == MomentId::mxxy)
        return (cx * cx - Stencil::cs2) * cy;
    else if constexpr (Id == MomentId::mxyy)
        return cx * (cy * cy - Stencil::cs2);
    else if constexpr (Id == MomentId::myyy)
        return cy * cy * cy - r::three * Stencil::cs2 * cy;
    else if constexpr (Id == MomentId::mxxyy)
    {
        return (cx * cx - Stencil::cs2) * (cy * cy - Stencil::cs2);
    }
    else
    {
        static_assert(always_false_v<Id>, "Unsupported MomentId in hermite().");
        return r::zero;
    }
}

template <auto Id>
__host__ __device__ __forceinline__
    real_t
    hermite(int i)
{
    return hermite<Id>(r_cast(Stencil::cx(i)), r_cast(Stencil::cy(i)));
}
