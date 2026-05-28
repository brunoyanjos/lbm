#pragma once

#include "../../core/types.cuh"
#include "../stencil_active.cuh"
#include "moment_id.cuh"

#include "lbm/meta/always_false.cuh"

template <auto Id>
struct ScaleFactor;

template <>
struct ScaleFactor<MomentId::rho>
{
    static constexpr real_t value = r::one;
};

template <>
struct ScaleFactor<MomentId::ux>
{
    static constexpr real_t value = Stencil::as2;
};

template <>
struct ScaleFactor<MomentId::uy>
{
    static constexpr real_t value = Stencil::as2;
};

template <>
struct ScaleFactor<MomentId::mxx>
{
    static constexpr real_t value = Stencil::as4 * r::half;
};

template <>
struct ScaleFactor<MomentId::mxy>
{
    static constexpr real_t value = Stencil::as4;
};

template <>
struct ScaleFactor<MomentId::myy>
{
    static constexpr real_t value = Stencil::as4 * r::half;
};

template <>
struct ScaleFactor<MomentId::mxxx>
{
    static constexpr real_t value = Stencil::as6 * r::sixth;
};

template <>
struct ScaleFactor<MomentId::mxxy>
{
    static constexpr real_t value = Stencil::as6 * r::half;
};

template <>
struct ScaleFactor<MomentId::mxyy>
{
    static constexpr real_t value = Stencil::as6 * r::half;
};

template <>
struct ScaleFactor<MomentId::myyy>
{
    static constexpr real_t value = Stencil::as6 * r::sixth;
};

template <auto Id>
__host__ __device__ __forceinline__ real_t
scale_factor()
{
    return ScaleFactor<Id>::value;
}

template <auto Id>
__host__ __device__ __forceinline__ real_t
inv_scale_factor()
{
    return r::one / scale_factor<Id>();
}
