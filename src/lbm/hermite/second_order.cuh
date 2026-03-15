#pragma once

#include "../../core/types.cuh"

namespace Hermite
{

    template <class StencilT>
    __host__ __device__ __forceinline__
        real_t
        Hxx(real_t cx)
    {
        return cx * cx - StencilT::cs2;
    }

    __host__ __device__ __forceinline__
        real_t
        Hxy(real_t cx, real_t cy)
    {
        return cx * cy;
    }

    template <class StencilT>
    __host__ __device__ __forceinline__
        real_t
        Hyy(real_t cy)
    {
        return cy * cy - StencilT::cs2;
    }

}