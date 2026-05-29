#pragma once

#include "lbm/meta/always_false.cuh"

#include "../../../core/types.cuh"
#include "../../stencil_active.cuh"
#include "../../moment/moment_id.cuh"
#include "../../hermite/hermite.cuh"

template <auto Id>
__host__ __device__ __forceinline__
    real_t
    moment_prefactor()
{
    if constexpr (Id == MomentId::rho)
    {
        return r::one;
    }
    else if constexpr (Id == MomentId::ux || Id == MomentId::uy)
    {
        return Stencil::as2;
    }
    else if constexpr (Id == MomentId::mxy)
    {
        return Stencil::as4;
    }
    else if constexpr (Id == MomentId::mxx || Id == MomentId::myy)
    {
        return r::half * Stencil::as4;
    }
    else if constexpr (Id == MomentId::mxxy || Id == MomentId::mxyy)
    {
        return r::half * Stencil::as6;
    }
    else if constexpr (Id == MomentId::mxxx || Id == MomentId::myyy)
    {
        return r::sixth * Stencil::as6;
    }
    else
    {
        static_assert(always_false_v<Id>, "Unsupported MomentId in moment_prefactor().");
        return r::zero;
    }
}

template <auto IdA, auto IdB>
__host__ __device__ __forceinline__
    real_t
    moment_factor(int i)
{
    return Stencil::w(i) * moment_prefactor<IdA>() * hermite<IdA>(i) * hermite<IdB>(i);
}