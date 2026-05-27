#pragma once
#include <cstdint>
#include "../stencil_active.cuh"
#include "../../core/indexing.cuh"

static_assert(Stencil::Q > 0, "Stencil::Q must be > 0");
static_assert(Stencil::Q <= 64, "mask_t supports only up to 64 stencil directions");

__host__ __device__ __forceinline__ static constexpr mask_t full_mask()
{
    static_assert(Stencil::Q <= 64, "Stencil::Q must be <= 64 for mask_t");

    if constexpr (Stencil::Q == sizeof(mask_t) * 8)
        return ~mask_t(0);
    else
        return (mask_t(1) << Stencil::Q) - mask_t(1);
}

__host__ __device__ __forceinline__ constexpr mask_t bit(int i)
{
    return (mask_t(1) << i);
}

constexpr __host__ __device__ __forceinline__ bool is_full_mask(mask_t m)
{
    const mask_t fm = full_mask();
    return (m & fm) == fm;
}
__host__ __device__ __forceinline__ bool dir_valid(mask_t valid_mask, int i)
{
    return (valid_mask & (mask_t(1) << i)) != mask_t(0);
}

__host__ __device__ __forceinline__ mask_t mask_opp(mask_t m)
{
    mask_t r = mask_t(0);

#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
    for (int i = 0; i < Stencil::Q; ++i)
    {
        if ((m & (mask_t(1) << i)) != mask_t(0))
            r |= (mask_t(1) << Stencil::opp(i));
    }

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

__host__ __device__ __forceinline__ int count_on_bits(mask_t m)
{
    m &= full_mask();

#if defined(__CUDA_ARCH__)
    if constexpr (sizeof(mask_t) == sizeof(uint32_t))
        return __popc(static_cast<uint32_t>(m));
    else
        return __popcll(static_cast<unsigned long long>(m));
#else
#if defined(__GNUC__) || defined(__clang__)
    if constexpr (sizeof(mask_t) == sizeof(uint32_t))
        return __builtin_popcount(static_cast<uint32_t>(m));
    else
        return __builtin_popcountll(static_cast<unsigned long long>(m));
#else
    int c = 0;
    while (m)
    {
        m &= (m - mask_t(1));
        ++c;
    }
    return c;
#endif
#endif
}

__host__ __device__ __forceinline__ int count_valid_dirs(mask_t valid_mask)
{
    return count_on_bits(valid_mask);
}

__host__ __device__ __forceinline__ int count_missing_dirs(mask_t valid_mask)
{
    return int(Stencil::Q) - count_on_bits(valid_mask);
}
