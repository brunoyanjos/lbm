#pragma once
#include <cstdint>
#include "../stencil_active.cuh"

__device__ __forceinline__ bool dir_valid(uint32_t valid_mask, int i)
{
    return (valid_mask & (uint32_t(1) << uint32_t(i))) != 0u;
}

__device__ __forceinline__ uint32_t mask_opp(uint32_t m)
{
    uint32_t r = 0u;
#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
        if (m & (1u << i))
            r |= (1u << Stencil::opp(i));
    return r;
}
