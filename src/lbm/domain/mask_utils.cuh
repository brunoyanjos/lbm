#pragma once
#include <cstdint>
#include "../stencil_active.cuh"
#include "../../core/indexing.cuh"
#include "../../core/types.cuh"
#include "domain_tags.cuh"

static_assert(Stencil::Q > 0, "Stencil::Q must be > 0");
static_assert(Stencil::Q <= 32, "valid_mask is uint32_t; use uint64_t if Q > 32");

constexpr __host__ __device__ __forceinline__ uint32_t full_mask()
{
    return (Stencil::Q == 32) ? 0xFFFFFFFFu : (uint32_t(1u) << uint32_t(Stencil::Q)) - 1u;
}

constexpr __host__ __device__ __forceinline__ bool is_full_mask(uint32_t m)
{
    const uint32_t fm = full_mask();
    return (m & fm) == fm;
}
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

__host__ __device__ __forceinline__
    uint8_t
    get_node_safe(const uint8_t *__restrict__ nodes, int x, int y)
{
    if ((unsigned)x >= (unsigned)NX || (unsigned)y >= (unsigned)NY)
        return to_u8(NodeId::SOLID);

    return nodes[idxGlobal(x, y)];
}

__host__ __device__ __forceinline__ int count_on_bits(uint32_t m)
{
    m &= full_mask();

#if defined(__CUDA_ARCH__)
    return __popc(m);
#else
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(m);
#else
    int c = 0;
    while (m)
    {
        m &= (m - 1u);
        ++c;
    }
    return c;
#endif
#endif
}

__host__ __device__ __forceinline__ int count_valid_dirs(uint32_t valid_mask)
{
    return count_on_bits(valid_mask);
}

__host__ __device__ __forceinline__ int count_missing_dirs(uint32_t valid_mask)
{
    return int(Stencil::Q) - count_on_bits(valid_mask);
}