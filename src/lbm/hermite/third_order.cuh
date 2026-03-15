#pragma once

#include "../../core/types.cuh"
#include "second_order.cuh"

namespace Hermite
{

    template <class StencilT>
    __host__ __device__ __forceinline__
        real_t
        Hxxy(real_t cx, real_t cy)
    {
        return Hxx<StencilT>(cx) * cy;
    }

    template <class StencilT>
    __host__ __device__ __forceinline__
        real_t
        Hxyy(real_t cx, real_t cy)
    {
        return cx * Hyy<StencilT>(cy);
    }

    template <class StencilT>
    __host__ __device__ __forceinline__
        real_t
        Hxxx(real_t cx)
    {
        return cx * (cx * cx - r::three * StencilT::cs2);
    }

    template <class StencilT>
    __host__ __device__ __forceinline__
        real_t
        Hyyy(real_t cy)
    {
        return cy * (cy * cy - r::three * StencilT::cs2);
    }

}